#!/usr/bin/env python3

from fuseq.option import Option
from fuseq.genomon import Genomon
from fuseq.blat import Blat
from fuseq.checker import Checker


def main():
    opt = Option()
    Checker.has_tools()
    genomon = Genomon(opt.refer())

    for mf_dir, mf_path in genomon.mf_dic.items():
        params = opt.copy()

        # Paths
        work_dir = f'{params.fuseq_root_dir}/{mf_dir}/_fuseq_work'
        fuseq_path = f'{params.fuseq_root_dir}/{mf_dir}/fusion_sequence.txt'

        # Add to params
        #params.mf_path = mf_path
        #params.jun_dic = genomon.jun_dic
        #params.work_dir = work_dir
        #params.fuseq_path = fuseq_path

        # Run Blat
        blat = Blat(mf_path, genomon.jun_dic, work_dir, fuseq_path, params)
        if params.is_restart:
            Checker.isdir(work_dir)
            blat.restart()
        else:
            blat.run()


if __name__ == '__main__':
    main()
