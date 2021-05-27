import subprocess
from fuseq.timer import Timer
from fuseq.base import Base

class Blat(Base):
    def __init__(self, params):
        super().__init__()
        self.params = params
        self.time_blat = Timer('Blat')

    def run(self):
        self.time_blat.start()

        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
blat {blat_opt} -noHead {reference} {inp_file} {out_file}
'''.format(work_dir=self.params.work_dir, blat_opt=self.params.blat_opt, reference=self.params.reference,
           inp_file=self.files['coll'], out_file=self.files['blat'])
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        _, err = p.communicate()
        if p.returncode != 0:
            print('[Error] at blat runtime')
            print(err)
            exit(1)

        self.time_blat.end()
        if self.params.print_time:
            self.time_blat.print()
