#!/usr/bin/env python3

from fuseq.option import Option
from fuseq.genomon import Genomon
from fuseq.blat import Blat
from fuseq.checker import Checker


def main():
    Checker.has_tools()
    opts = Option().create()
    genomon = Genomon(opts)
    for mf_dir, mf_path in genomon.mf_dic.items():
        work_dir = f'{opts.fuseq_root_dir}/{mf_dir}/_fuseq_work'
        fuseq_path = f'{opts.fuseq_root_dir}/{mf_dir}/fusion_sequence.txt'
        blat = Blat(mf_path, genomon.jun_dic, work_dir, fuseq_path, opts)
        if opts.restart_filter or opts.restart_blat:
            Checker.isdir(work_dir)
            if opts.restart_filter:
                blat.restart_from_filter()
            else:
                blat.restart_from_blat()
        else:
            blat.run()


if __name__ == '__main__':
    main()
