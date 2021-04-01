import argparse
import os
from fuseq import __version__
from fuseq.checker import Checker


class Option:
    def __init__(self):
        self._parse()
        self._create()

    def _parse(self):
        prog = 'fuseq'
        parser = argparse.ArgumentParser(prog=prog)
        parser.add_argument('genomon_dir', metavar='input_directory', action='store', help='input directory')
        parser.add_argument('--version', action='version', version=f'{prog}: {__version__}')
        self._args = parser.parse_args()

    def _create(self):
        args = self._args
        Checker.isdir(args.genomon_dir)

        self._opts = {}
        self._opts['genomon_dir'] = os.path.abspath(args.genomon_dir)

    @property
    def opts(self):
        return self._opts
