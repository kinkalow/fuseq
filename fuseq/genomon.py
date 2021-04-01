import glob
import os
from fuseq.checker import Checker

class Genomon:
    def __init__(self, opts):
        self.genomon_dir = opts.genomon_dir
        self.use_filt = opts.use_filt

        self._mf_path = {}
        self._jun_path = {}
        self._paths = [self._mf_path, self._jun_path]

        # Create required file paths and check for their existence
        self._create_mf_path()
        self._create_jun_path()

    def _create_mf_path(self):
        f = 'merge_fusionfusion_filt.txt' if self.use_filt else 'merge_fusionfusion.txt'
        base = f'{self.genomon_dir}/post_analysis'
        Checker.isdir(base)
        for d in os.listdir(base):
            path = f'{base}/{d}/{f}'
            Checker.isfile(path)
            self._mf_path[d] = path

    def _create_jun_path(self):
        ext1 = '.junction'
        ext2 = '.sam'
        base = f'{self.genomon_dir}/star'
        Checker.isdir(base)
        for d in os.listdir(base):
            li = glob.glob(f'{base}/{d}/*{ext1}')
            Checker.isonefile(li)
            f_jun = li[0]
            f_sam = f_jun[:-len(ext1)] + ext2
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
