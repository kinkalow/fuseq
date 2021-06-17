#!/usr/bin/env python3

from fuseq.option import Option
from fuseq.genomon import Genomon
from fuseq.pipeline import Pipeline


def main():

    opt = Option()
    genomon = Genomon(opt.refer())

    for mf_dir, mf_path in genomon.mf_dic.items():
        params = opt.copy()

        # Paths
        work_dir = f'{params.fuseq_root_dir}/{mf_dir}/{params.work_dirname}'
        swork_dir = f'{params.fuseq_root_dir}/{mf_dir}/{params.work_dirname}/{params.swork_dirname}'
        fuseq_path = f'{params.fuseq_root_dir}/{mf_dir}/{params.fuseq_filename}'
        inputs = {'mf_path': mf_path, 'star_dir': genomon.star_dir}

        # Add to params
        params.work_dir = work_dir
        params.swork_dir = swork_dir
        params.fuseq_path = fuseq_path
        params.inputs = inputs

        # Run
        pipeline = Pipeline(params)
        pipeline.run()


if __name__ == '__main__':
    main()
