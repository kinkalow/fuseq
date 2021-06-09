import re
import subprocess
import time

class Base:
    def __init__(self):
        self.files = {'params': 'params', 'coll': 'collect', 'blat': 'blat',
                      'filtmatch': 'filter_match', 'filtmiss': 'filter_miss',
                      'filtwar': 'filter_warning', 'filterr': 'filter_error',
                      'coll_res': 'collect_restart', 'blat_res': 'blat_restart'}

    def _run_cmd(self, cmd, name=None, ignore_err=False):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        out = out.decode().rstrip('\n')
        err = err.decode().rstrip('\n')
        if ignore_err:
            return out, err, p.returncode
        if p.returncode != 0:
            if not name:
                print(f'[Error] at {name} runtime')
            print(err)
            exit(1)
        return out

    # Array job
    def _run_jobs(self, code, path, num_parallels, name=None):
        '''Execute array job and wait for completion
           Check the return code'''

        with open(path, 'w') as f:
            f.write(code)

        cmd = f'qsub -terse -sync y -t 1-{num_parallels}:1 {path}'
        out, err, ret = self._run_cmd(cmd, 'qsub', ignore_err=True)
        # out contains jobid and return codes

        if err:
            print('[stderr]')
            print(err)
            print('[stdout]')
            print(out)
            exit(1)

        # Obtain return codes from qsub
        num_matched = 0
        is_okay = True
        for o in out.split('\n')[1:]:
            match = re.match(r'^Job (.*) exited with exit code ([\d]+)\.$', o)
            if match:
                num_matched += 1
                jobid = match.group(1)
                ret_code = match.group(2)
                if ret_code != '0':
                    print(f'[Error] {jobid} failed with exit code {ret_code}')
                    is_okay = False
        if not is_okay:
            if num_matched != num_parallels:
                print('[Error] Cannot get all exit codes')
            exit(1)

        # Obtain return codes from qacct
        if num_matched != num_parallels:
            print('[Warning] Output message changed after job execution')
            jobid = out.split('\n')[0]
            time.sleep(10)
            for _ in range(20):
                cmd = f'qacct -j {jobid} | grep exit_status'
                out, _, ret = self._run_cmd(cmd, 'qacct', ignore_err=True)
                if ret == 0:
                    break
                else:
                    time.sleep(30)
            else:
                print('[Error] Cannot get exit_status using qacct')
                exit(1)
            try:
                rets = [int(re.match(r'^exit_status.*([\d]+)', o).group(1)) for o in out.split('\n')]
                if len(rets) != num_parallels:
                    raise Exception
            except Exception:
                print(f'[Error] Cannot get exit status of jobid {jobid}')
                exit(1)
            for ret in rets:
                if ret != 0:
                    print(f'[Error] Job {jobid} failed with exit codes str(rets)')
                    exit(1)
