import glob
import os
from fuseq.checker import Checker

class Genomon:
    def __init__(self, genomon_dir):
        self.genomon_dir = genomon_dir

        self._mf_path = {}
        self._jun_path = {}
        self._paths = [self._mf_path, self._jun_path]

        self._create_mf_path()
        self._create_jun_path()

    def _create_mf_path(self):
        f = 'merge_fusionfusion.txt'
        #f = 'merge_fusionfusion_filt.txt'
        base = f'{self.genomon_dir}/post_analysis'
        Checker.isdir(base)
        for d in os.listdir(base):
            path = f'{base}/{d}/{f}'
            Checker.isfile(path)
            self._mf_path[d] = path

    def _create_jun_path(self):
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
            self._jun_path[d] = f_jun

    @property
    def paths(self):
        return self._paths

    @property
    def mf_dic(self):
        return self._mf_path

    @property
    def jun_dic(self):
        return self._jun_path
