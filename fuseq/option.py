import argparse
import os
from fuseq import __version__
from fuseq.checker import Checker


class Option:
    def __init__(self):
        self._parse()
        self._check()

    def _parse(self):
        prog = 'fuseq'
        parser = argparse.ArgumentParser(prog=prog)
        parser.add_argument('genomon_dir', metavar='input_directory', type=str, help='input directory')
        parser.add_argument('output_dir', metavar='output_directory', type=str, help='output directory')
        parser.add_argument('--no-use-filt', default=False, action='store_true', help='use merge_fusionfusion.txt instead of merge_fusionfusion_filt.txt')
        parser.add_argument('--reference', default='/share/pub/genomon/.genomon_local/genomon_pipeline-2.6.3/database/GRCh37/GRCh37.fa', type=str, help='reference path')
        parser.add_argument('--version', action='version', version=f'{prog}: {__version__}')
        self.args = parser.parse_args()

    def _check(self):
        args = self.args
        Checker.isdir(args.genomon_dir)
        Checker.isfile(args.reference)
        #if os.path.isdir(args.output_dir):
        #    print('output directory exists')
        #    exit(1)

    def create(self):
        args = self.args
        opts = {}
        opts['genomon_dir'] = os.path.abspath(args.genomon_dir)
        opts['out_base_dir'] = os.path.abspath(args.output_dir)
        opts['preprocess_parallel_num'] = 4
        opts['reference'] = os.path.abspath(args.reference)
        opts['use_filt'] = False if args.no_use_filt else True
        opts = argparse.Namespace(**opts)
        return opts
