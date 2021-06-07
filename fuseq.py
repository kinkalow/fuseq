#!/usr/bin/env python3

from fuseq.option import Option
from fuseq.genomon import Genomon
from fuseq.pipeline import Pipeline
from fuseq.checker import Checker


def main():

    opt = Option()
    Checker.has_tools()
    genomon = Genomon(opt.refer())

    for mf_dir, mf_path in genomon.mf_dic.items():
        params = opt.copy()

        # Paths
        work_dir = f'{params.fuseq_root_dir}/{mf_dir}/{params.work_dirname}'
        skwork_dir = f'{params.fuseq_root_dir}/{mf_dir}/{params.work_dirname}/{params.skwork_dirname}'
        fuseq_path = f'{params.fuseq_root_dir}/{mf_dir}/{params.fuseq_filename}'
        inputs = {'mf_path': mf_path, 'star_dir': genomon.star_dir}

        # Add to params
        params.work_dir = work_dir
        params.skwork_dir = skwork_dir
        params.fuseq_path = fuseq_path
        params.inputs = inputs
        params.num_blat_parallels = 100  # FIXME:

        # Run Blat
        pipeline = Pipeline(params)
        if params.is_restart:
            Checker.isdir(work_dir)
            pipeline.restart()
        elif params.is_shirokane:
            pipeline.run_on_shirokane()
        else:
            pipeline.run()


if __name__ == '__main__':
    main()
