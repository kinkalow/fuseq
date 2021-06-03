import glob
from fuseq.timer import Timer
from fuseq.base import Base

class Blat(Base):
    def __init__(self, params):
        super().__init__()
        self.params = params
        self.parallel_number = self.params.num_parallel_blat

    #def __create_blat_cmd(self, work_dir, inp_file, out_file):
        #cmd = '''\
#!/bin/bash
#set -eu
#cd {work_dir}
#blat {blat_opt} -noHead {reference} {inp_file} {out_file}
#'''.format(work_dir=work_dir, inp_file=inp_file, out_file=out_file,
        #blat_opt=self.params.blat_opt, reference=self.params.reference)
        #return cmd

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
    # On shirokane
    #

    def __split(self):
        prefix = self.files['coll']
        cmd = '''\
#!/bin/bash
set -eu

cd {work_dir}
line_cnt=$(wc -l {inp_file} | cut -f 1 -d ' ')
max_proc=$((line_cnt/2))
number=$([ {parallel_number} -lt $max_proc ] && echo {parallel_number} || echo $max_proc)
length=${{#number}}
lines=$([ $((line_cnt/2%number)) -eq 0 ] && echo $((line_cnt/number)) || echo $(((line_cnt/2/number+1)*2)))

cd {skwork_dir}
split -a $length -d -l $lines --numeric-suffixes=1 ../{inp_file} {prefix}
'''.format(work_dir=self.params.work_dir, skwork_dir=self.params.skwork_dir,
           inp_file=self.files['coll'], parallel_number=self.parallel_number, prefix=prefix)
        self._run_cmd(cmd, 'split')
        # Update the parallel number of blat
        self.parallel_number = len(glob.glob(f'{self.params.skwork_dir}/{self.files["coll"]}*'))

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
           length=len(list(str(self.parallel_number))),
           blat_path=blat_path, blat_opt=self.params.blat_opt,
           reference=self.params.reference, inp_file=self.files['coll'])
        out_file = f'{self.params.skwork_dir}/{self.files["blat"]}.sh'
        with open(out_file, 'w') as f:
            f.write(script)
        self._run_cmd_on_shirokane(out_file, self.parallel_number, 'blat_shirokane')

    def __concat(self):
        length = len(list(str(self.parallel_number)))
        cmd = '''\
!/bin/bash
set -eu
cd {skwork_dir}
cat {inp_files} > ../{out_file}
'''.format(skwork_dir=self.params.skwork_dir,
           inp_files=' '.join([f'{self.files["blat"]}{str(i).zfill(length)}'
                               for i in range(1, self.parallel_number + 1)]),
           out_file=self.files['blat'])
        self._run_cmd(cmd, 'cat_blat')

    @Timer('blat')
    def run_batch(self):
        self.__split()
        self.__blat()
        self.__concat()
