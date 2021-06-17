import glob
import os
from fuseq.checker import Checker

class Genomon:
    def __init__(self, opts):
        self.opts = opts
        # fusion path and star directory path
        self.__mf_path = {}
        self.__star_dir = opts.star_dir if opts.star_dir else f'{opts.genomon_root_dir}/star'
        # Create fusion paths and check if files exist
        self.__create_mf_path()
        self.__check_inside_star()

    def __create_mf_path(self):
        if self.opts.fusion_file:
            path = self.opts.fusion_file
            Checker.isfile(path)
            self.__mf_path['noname'] = path
        else:
            f = 'merge_fusionfusion_filt.txt' if self.opts.use_filt else 'merge_fusionfusion.txt'
            base = f'{self.opts.genomon_root_dir}/post_analysis'
            Checker.isdir(base)
            for d in os.listdir(base):
                path = f'{base}/{d}/{f}'
                Checker.isfile(path)
                self.__mf_path[d] = path

    def __check_inside_star(self):
        # Check if star direcotry exits
        star_dir = self.__star_dir
        Checker.isdir(star_dir)
        # Check for the existence of junction and sam files
        for d in os.listdir(star_dir):
            li = glob.glob(f'{star_dir}/{d}/*.junction')
            Checker.isonefile(li)
            f_jun = li[0]
            f_sam = f_jun[:-8] + 'sam'
            Checker.isfile(f_sam)

    @property
    def mf_dic(self):
        return self.__mf_path

    @property
    def star_dir(self):
        return self.__star_dir
