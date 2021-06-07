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
        self.num_numeric_suffix = len(list(str(self.num_parallels)))

    def __calculate_numbers(self):
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
num_lines=$(wc -l {inp_file} | cut -f 1 -d ' ')
echo -n $num_lines
'''.format(work_dir=self.params.work_dir, inp_file=self.files['coll'])
        num_lines = int(self._run_cmd(cmd, 'num_lines'))
        maxnum_parallels = int(num_lines / 2)
        num_parallels = self.params.num_blat_parallels if self.params.num_blat_parallels < maxnum_parallels else maxnum_parallels
        return num_parallels, num_lines

    def __split(self):
        ttl_read = int(self.num_coll_lines / 2)
        n_read_per_file = int(ttl_read // self.num_parallels)
        n_plus1_file = ttl_read % self.num_parallels
        n_plus0_file = self.num_parallels - n_plus1_file
        n_line1_per_file = 2 * (n_read_per_file + 1)
        n_line0_per_file = 2 * n_read_per_file
        prefix = self.files['coll']

        if n_plus1_file == 0:
            cmd = '''\
#!/bin/bash
set -eu
cd {skwork_dir}
split -a {length} -d -l {lines} --numeric-suffixes=1 ../{inp_file} {prefix}
'''.format(skwork_dir=self.params.skwork_dir, length=self.num_numeric_suffix,
           lines=n_line1_per_file, inp_file=self.files['coll'], prefix=prefix)
            self._run_cmd(cmd, 'split_coll1')
        else:
            num_plus1_lines = n_line1_per_file * n_plus1_file
            num_plus0_lines = n_line0_per_file * n_plus0_file
            cmd = '''\
#!/bin/bash
set -eu
cd {skwork_dir}
head -{num_plus1_lines} ../{inp_file} | split -a {length} -d -l {lines1} --numeric-suffixes=1 - {prefix}
tail -{num_plus0_lines} ../{inp_file} | split -a {length} -d -l {lines0} --numeric-suffixes={n_suf0} - {prefix}
'''.format(skwork_dir=self.params.skwork_dir, num_plus1_lines=num_plus1_lines,
           inp_file=self.files['coll'], length=self.num_numeric_suffix,
           lines1=n_line1_per_file, prefix=prefix, num_plus0_lines=num_plus0_lines,
           lines0=n_line0_per_file, n_suf0=n_plus1_file + 1)
            self._run_cmd(cmd, 'split_coll2')

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
num=$(printf "%0{length}d" ${{SGE_TASK_ID}})
cd {skwork_dir}
{blat_path} {blat_opt} -noHead {reference} {inp_file}${{num}} {out_file}${{num}}
'''.format(skwork_dir=self.params.skwork_dir, out_file=self.files['blat'],
           length=self.num_numeric_suffix,
           blat_path=blat_path, blat_opt=self.params.blat_opt,
           reference=self.params.reference, inp_file=self.files['coll'])
        out_file = f'{self.params.skwork_dir}/{self.files["blat"]}.sh'
        with open(out_file, 'w') as f:
            f.write(script)
        self._run_cmd_on_shirokane(out_file, self.num_parallels, 'blat_shirokane')

    def __concat(self):
        length = self.num_numeric_suffix
        cmd = '''\
!/bin/bash
set -eu
cd {skwork_dir}
cat {inp_files} > ../{out_file}
'''.format(skwork_dir=self.params.skwork_dir,
           inp_files=' '.join([f'{self.files["blat"]}{str(i).zfill(length)}'
                               for i in range(1, self.num_parallels + 1)]),
           out_file=self.files['blat'])
        self._run_cmd(cmd, 'cat_blat')

    @Timer('blat')
    def run_batch(self):
        self.__split()
        self.__blat()
        self.__concat()
