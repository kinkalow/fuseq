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
        parser.add_argument('genomon_dir', metavar='input_directory', action='store', help='input directory')
        parser.add_argument('output_dir', metavar='output_directory', help='output directory')
        parser.add_argument('--version', action='version', version=f'{prog}: {__version__}')
        self.args = parser.parse_args()

    def _check(self):
        args = self.args
        Checker.isdir(args.genomon_dir)
        #if os.path.isdir(args.output_dir):
        #    print('output directory exists')
        #    exit(1)

    def create(self):
        args = self.args
        opts = {}
        opts['genomon_dir'] = os.path.abspath(args.genomon_dir)
        opts['out_dir'] = os.path.abspath(args.output_dir)
        opts = argparse.Namespace(**opts)
        return opts
