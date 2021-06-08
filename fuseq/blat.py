from fuseq.timer import Timer
from fuseq.base import Base

class Blat(Base):
    def __init__(self, params):
        super().__init__()
        self.params = params

    @Timer('blat')
    def run(self):
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
blat {blat_opt} -noHead {reference} {inp_file} {out_file}
'''.format(work_dir=self.params.work_dir, blat_opt=self.params.blat_opt,
           reference=self.params.reference, inp_file=self.files['coll'],
           out_file=self.files['blat'])
        self._run_cmd(cmd, 'blat')


#
# Parallelization of Blat on Shirokane
#

class PBlat(Base):
    def __init__(self, params):
        super().__init__()
        self.params = params
        self.num_parallels, self.num_coll_lines = self.__calculate_numbers()
        self.num_numeric_suffixes = len(list(str(self.num_parallels)))

    def __calculate_numbers(self):
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
n_lines=$(wc -l {inp_file} | cut -f 1 -d ' ')
echo -n $n_lines
'''.format(work_dir=self.params.work_dir, inp_file=self.files['coll'])
        n_lines = int(self._run_cmd(cmd, 'n_lines'))
        max_parallels = int(n_lines / 2)
        n_parallels = self.params.num_blat_parallels if self.params.num_blat_parallels < max_parallels else max_parallels
        return n_parallels, n_lines

    def __split(self):
        n_reads = int(self.num_coll_lines / 2)
        n_reads0 = int(n_reads // self.num_parallels)
        n_files1 = n_reads % self.num_parallels   # Number of files with lines equal to read number plus 1
        n_files0 = self.num_parallels - n_files1  # Number of files with lines equal to read number plus 0
        n_lines1_per_file = 2 * (n_reads0 + 1)
        n_lines0_per_file = 2 * n_reads0
        n_lines1 = n_lines1_per_file * n_files1
        n_lines0 = n_lines0_per_file * n_files0
        prefix = self.files['coll']

        cmd = '''\
#!/bin/bash
set -eu
cd {skwork_dir}
head -{n_lines1} ../{inp_file} | split -a {length} -d -l {lines1} --numeric-suffixes=1 - {prefix}
tail -{n_lines0} ../{inp_file} | split -a {length} -d -l {lines0} --numeric-suffixes={numr_sfx} - {prefix}
'''.format(skwork_dir=self.params.skwork_dir, n_lines1=n_lines1,
           inp_file=self.files['coll'], length=self.num_numeric_suffixes,
           lines1=n_lines1_per_file, prefix=prefix, n_lines0=n_lines0,
           lines0=n_lines0_per_file, numr_sfx=n_files1 + 1)
        self._run_cmd(cmd, 'split_coll')

    def __blat(self):
        blat_path = self._run_cmd('which blat', 'which_blat')
        script = '''\
#!/usr/local/bin/nosh
#$ -S /usr/local/bin/nosh
#$ -cwd
#$ -l s_vmem=8G,mem_req=8G
#$ -e {skwork_dir}/{out_file}.log
#$ -o {skwork_dir}/{out_file}.log
set -eu
id=$(printf "%0{length}d" ${{SGE_TASK_ID}})
cd {skwork_dir}
{blat_path} {blat_opt} -noHead {reference} {inp_file}${{id}} {out_file}${{id}}
'''.format(skwork_dir=self.params.skwork_dir, out_file=self.files['blat'],
           length=self.num_numeric_suffixes,
           blat_path=blat_path, blat_opt=self.params.blat_opt,
           reference=self.params.reference, inp_file=self.files['coll'])
        out_file = f'{self.params.skwork_dir}/{self.files["blat"]}.sh'
        with open(out_file, 'w') as f:
            f.write(script)
        self._run_cmd_on_shirokane(out_file, self.num_parallels, 'blat_shirokane')

    def __concat(self):
        cmd = '''\
!/bin/bash
set -eu
cd {skwork_dir}
cat {inp_files} > ../{out_file}
'''.format(skwork_dir=self.params.skwork_dir,
           inp_files=' '.join([f'{self.files["blat"]}{str(i).zfill(self.num_numeric_suffixes)}'
                               for i in range(1, self.num_parallels + 1)]),
           out_file=self.files['blat'])
        self._run_cmd(cmd, 'cat_blat')

    @Timer('blat')
    def run_batch(self):
        self.__split()
        self.__blat()
        self.__concat()
