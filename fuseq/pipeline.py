import json
import os
import shutil
from fuseq.checker import Checker
from fuseq.base import Base
from fuseq.collection import Collection
from fuseq.blat import Blat, PBlat
from fuseq.blat_filter import BlatFilter

class Pipeline(Base):

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.breakinfo_path = f'{self.params.work_dir}/{self.files["breakinfo"]}'

    #
    # Utility
    #

    def __save_params(self):
        path = f'{self.params.work_dir}/{self.files["params"]}'
        if os.path.exists(path):
            base_path = path
            num = 1
            while True:
                path = f'{base_path}.{num}'
                if not os.path.exists(path):
                    break
                num += 1
        with open(path, 'w') as f:
            d = vars(self.params)
            max_len = max([len(k) for k in d.keys()])
            for key, val in sorted(d.items(), key=lambda x: x[0]):
                space = ' ' * (max_len - len(key))
                f.write(f'{key}{space}: {val}\n')

    def __save_breakinfo(self, breakinfo):
        with open(self.breakinfo_path, 'w') as f:
            json.dump(breakinfo, f)

    def __load_breakinfo(self):
        Checker.isfile(self.breakinfo_path)
        with open(self.breakinfo_path) as f:
            breakinfo = json.load(f)
        return breakinfo

    def __delete_work_dir(self, make_empty_dir=False):
        shutil.rmtree(self.params.work_dir, ignore_errors=True)
        if make_empty_dir:
            os.makedirs(self.params.work_dir)
            os.makedirs(self.params.swork_dir)

    #
    # Restart
    #

    def __on_blat_restart(self, breakinfo):
        coll_inp = f"{self.params.work_dir}/{self.files['coll']}"
        coll_out = f"{self.params.work_dir}/{self.files['coll_res']}"

        # Conditions for obtaining targets
        if self.params.readname_filt:
            target = self.params.readname_filt + '\n'

            def cond():
                return readname[readname.index('_') + 1:] == target
        else:
            target = self.params.seq_filt + '\n'

            def cond():
                return seq == target

        # Filter collection data
        linenrs = []
        cnts = {}
        with open(coll_out, 'w') as fw:
            with open(coll_inp, 'r') as fr:
                while True:
                    readname = fr.readline()
                    if not readname:
                        break
                    seq = fr.readline()
                    if cond():
                        linenr = int(readname[1:readname.index('_')].split('-')[0])  # >2-1_READNAME => 2
                        linenrs.append(linenr)
                        cnts[linenr] = cnts[linenr] + 1 if linenr in cnts else 1
                        fw.write(readname)
                        fw.write(seq)

        # Check
        if not linenrs:
            print('[Error] No input data after filtering')
            exit(1)

        # Filter breakinfo
        linenrs = set(linenrs)
        breakinfo = [b for b in breakinfo if int(b['linenr']) in linenrs]
        for b in breakinfo:
            b['cnt'] = cnts[b['linenr']]

        return breakinfo

    def __on_blat_filter_restart(self, breakinfo):
        breakinfo = self.__on_blat_restart(breakinfo)

        # Get readnames
        coll_inp = f"{self.params.work_dir}/{self.files['coll_res']}"
        readnames = []
        with open(coll_inp, 'r') as f:
            while True:
                readname = f.readline().rstrip('\n')
                if not readname:
                    break
                readnames.append(readname[1:])
                f.readline()

        # Filter blat data
        blat_inp = f"{self.params.work_dir}/{self.files['blat']}"
        blat_out = f"{self.params.work_dir}/{self.files['blat_res']}"
        with open(blat_inp, 'r') as fr:
            with open(blat_out, 'w') as fw:
                for blat in fr:
                    if blat.split('\t')[9] in readnames:
                        fw.write(blat)

        return breakinfo

    def __filter_on_restart(self, breakinfo):
        if not self.params.readname_filt and not self.params.seq_filt:
            return breakinfo

        if self.params.restart_blat:
            breakinfo = self.__on_blat_restart(breakinfo)
        else:
            breakinfo = self.__on_blat_filter_restart(breakinfo)

        # Change file names
        self.files['coll'] = self.files['coll_res']
        self.files['blat'] = self.files['blat_res']

        return breakinfo

    #
    # Start
    #

    def run(self, restart=False):
        if self.params.is_restart:
            # Preprocess
            Checker.isdir(self.params.work_dir)
            self.__save_params()
            breakinfo = self.__load_breakinfo()
            breakinfo = self.__filter_on_restart(breakinfo)

        else:
            # Preprocess
            self.__delete_work_dir(make_empty_dir=True)
            self.__save_params()
            # Collection
            breakinfo = Collection(self.params).run()
            self.__save_breakinfo(breakinfo)

        if self.params.stop_blat:
            return

        # Blat
        if not self.params.restart_filter:
            if self.params.on_shirokane:
                PBlat(self.params).run()
            else:
                Blat(self.params).run()

        if self.params.stop_filter:
            return

        # Filter
        BlatFilter(self.params, breakinfo).run()

        # Postprocess
        if self.params.delete_work:
            self.__delete_work_dir()
