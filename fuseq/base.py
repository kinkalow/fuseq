import subprocess

class Base:
    def __init__(self):
        self.files = {'params': 'params', 'coll': 'collect', 'blat': 'blat',
                      'filtmatch': 'filter_match', 'filtmiss': 'filter_miss',
                      'coll_res': 'collect_restart', 'blat_res': 'blat_restart'}

    def run_cmd(self, cmd, name=None, ignore_err=False):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        if ignore_err:
            return out.decode(), err.decode(), p.returncode
        if p.returncode != 0:
            if not name:
                print('[Error] at blat runtime')
            print(err.decode())
            exit(1)
        return out.decode()
