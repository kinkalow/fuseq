import argparse
import os
from fuseq import __version__
from fuseq.checker import Checker


class Option:
    def __init__(self):
        self._parse()
        self._check()
        self._create()

    def _parse(self):
        prog = 'fuseq'

        #class DelWorkAct(argparse.Action):
        #    def __call__(self, parser, namespace, values, option_string=None):
        #        if self.dest == 'delete_work':
        #            delete = True
        #        elif self.dest == 'no_delete_work':
        #            delete = False
        #        if 'tmp_delete_work' not in namespace:
        #            setattr(namespace, 'tmp_delete_work', delete)
        #        else:
        #            namespace.tmp_delete_work = delete

        parser = argparse.ArgumentParser(prog=prog)
        parser.add_argument('genomon_root_dir', metavar='genomon_root_directory', type=str, help='Input directory')
        parser.add_argument('fuseq_root_dir', metavar='fuseq_root_directory', type=str, help='Output directory')
        parser.add_argument('-b', '--blat-opt', default='', type=str, help='Blat options')
        #parser.add_argument('-b', '--blat-opt', default='-minScore=20 -minMatch=1', type=str, help='Blat options')
        parser.add_argument('-d', '--star-dir', default='', type=str, help='Alternative star directory in genomon output')
        parser.add_argument('-f', '--fusion-file', default='', type=str, help='Alternative merge_fusionfusion_filt.txt file in genomon output')
        parser.add_argument('-p', '--coll-procs', default=4, type=int, help='Number of parallel processes for collecting data for Blat input')
        parser.add_argument('-w', '--no-delete-work', default=False, action='store_true', help='Do not delete work directory')
        #parser.add_argument('-w', '--delete-work', default=False, action='store_true', help='Delete work directory')
        #parser.add_argument('-w', '--delete-work', default=False, nargs=0, action=DelWorkAct, help='Delete work directory on restart')
        parser.add_argument('-B', '--restart-blat', default=False, action='store_true', help='Restart from Blat')
        parser.add_argument('-F', '--restart-filter', default=False, action='store_true', help='Restart from Blat filter')
        parser.add_argument('--no-use-filt', default=False, action='store_true', help='Use merge_fusionfusion.txt instead of merge_fusionfusion_filt.txt')
        #parser.add_argument('--no-delete-work', default=True, nargs=0, action=DelWorkAct, help='Do not delete work directory when not restarting')
        parser.add_argument('--reference', default='/share/pub/genomon/.genomon_local/genomon_pipeline-2.6.3/database/GRCh37/GRCh37.fa', type=str, help='Reference path')
        parser.add_argument('--end', default=1, type=int, help='Extend the end position of breakpoints at Blat filtering')
        parser.add_argument('--start', default=0, type=int, help='Extend the start position of breakpoints at Blat filtering')
        #parser.add_argument('--show-time', default=False, action='store_true', help='Display elapsed time')
        parser.add_argument('--version', action='version', version=f'{prog}: {__version__}')
        self.args = parser.parse_args()

    def _check(self):
        args = self.args
        Checker.isdir(args.genomon_root_dir)
        Checker.isfile(args.reference)

    def _create(self):
        args = self.args

        # Change Args
        args.genomon_root_dir = os.path.abspath(args.genomon_root_dir)
        args.fuseq_root_dir = os.path.abspath(args.fuseq_root_dir)
        args.fusion_file = os.path.abspath(args.fusion_file) if args.fusion_file else ''
        args.star_dir = os.path.abspath(args.star_dir) if args.star_dir else ''
        args.reference = os.path.abspath(args.reference)
        #
        args.delete_work = not args.no_delete_work
        args.use_filt = False if args.no_use_filt else True
        args.bp_start_extn = args.start
        args.bp_end_extn = args.end
        del args.no_delete_work
        del args.no_use_filt
        del args.start
        del args.end
        #
        # Change Arg
        #if 'tmp_delete_work' in args:
        #    delete_work = args.tmp_delete_work
        #    del args.tmp_delete_work
        #else:
        #    delete_work = args.delete_work if args.is_restart else args.no_delete_work
        #args.delete_work = delete_work

        # New Arg
        args.is_restart = True if args.restart_blat or args.restart_filter else False

    def refer(self):
        return self.args

    def copy(self):
        return argparse.Namespace(**vars(self.args))  # Deep copy
