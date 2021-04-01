import glob
import os
from fuseq.checker import Checker

class Genomon:
    def __init__(self, genomon_dir):
        self.genomon_dir = genomon_dir

        self._path_mf = {}
        self._path_jun = {}
        self._paths = [self._path_mf, self._path_jun]

        self._create_path_mf()
        self._create_path_jun()

    def _create_path_mf(self):
        f = 'merge_fusionfusion.txt'
        #f = 'merge_fusionfusion_filt.txt'
        base = f'{self.genomon_dir}/post_analysis'
        Checker.isdir(base)
        for d in os.listdir(base):
            path = f'{base}/{d}/{f}'
            Checker.isfile(path)
            self._path_mf[d] = path

    def _create_path_jun(self):
        f1 = '.junction'
        f2 = '.sam'
        base = f'{self.genomon_dir}/star'
        Checker.isdir(base)
        for d in os.listdir(base):
            li = glob.glob(f'{base}/{d}/*{f1}')
            Checker.isonefile(li)
            f_jun = li[0]
            f_sam = f_jun[:-len(f1)] + f2
            Checker.isfile(f_sam)
            self._path_jun[d] = f_jun

    @property
    def paths(self):
        return self._paths

    @property
    def path_mf(self):
        return self._path_mf

    @property
    def path_jun(self):
        return self._path_jun
