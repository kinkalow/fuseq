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
'''.format(work_dir=self.params.work_dir, blat_opt=self.params.blat_opt, reference=self.params.reference,
           inp_file=self.files['coll'], out_file=self.files['blat'])
        self.run_cmd(cmd, 'blat')
