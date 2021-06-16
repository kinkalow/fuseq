import argparse
import configparser
import os
import re
from fuseq import __version__
from fuseq.checker import Checker


class Option:
    def __init__(self):
        self.__parse()
        self.__check()
        self.__create()

    def __parse(self):
        prog = 'fuseq'

        parser = argparse.ArgumentParser(prog=prog)
        parser.add_argument('genomon_root_directory', type=str, help='Input directory')
        parser.add_argument('fuseq_root_directory', type=str, help='Output directory')
        parser.add_argument('-b', '--blat-options', default='', type=str, help='Blat options')
        parser.add_argument('-d', '--star-directory', default='', type=str, help='Alternative star directory in genomon output')
        parser.add_argument('-l', '--lines', default='', type=str, help='Line number in a fusion gene detection file')
        parser.add_argument('-f', '--fusion-file', default='', type=str, help='Alternative merge_fusionfusion_filt.txt file in genomon output')
        parser.add_argument('-w', '--no-delete-work', default=False, action='store_true', help='Do not delete work directory')
        parser.add_argument('-B', '--restart-blat', default=False, action='store_true', help='Restart from Blat')
        parser.add_argument('-F', '--restart-filter', default=False, action='store_true', help='Restart from Blat filter')
        #
        parser.add_argument('--collection-processes', default=4, type=int, help='Number of parallel processes for collection computation')
        parser.add_argument('--no-use-filt', default=False, action='store_true', help='Use merge_fusionfusion.txt instead of merge_fusionfusion_filt.txt')
        parser.add_argument('--readname', default='', type=str, help='Filtering with readname')
        parser.add_argument('--sequence', default='', type=str, help='Filtering with sequence')
        parser.add_argument('--start', default=0, type=int, help='Extend the start position of breakpoints at Blat filtering')
        parser.add_argument('--end', default=0, type=int, help='Extend the end position of breakpoints at Blat filtering')
        parser.add_argument('--reference', default='/share/pub/genomon/.genomon_local/genomon_pipeline-2.6.3/database/GRCh37/GRCh37.fa', type=str, help='Reference path')
        # Options on Shirokane
        parser.add_argument('--shirokane', default=False, action='store_true', help='Compute on Shirokane')
        parser.add_argument('--collection-tasks', default=100, type=int, help='Number of array job tasks for collection computation')
        parser.add_argument('--blat-tasks', default=100, type=int, help='Number of array job tasks for blat computation')
        #
        parser.add_argument('--no-check-position-interval', default=False, action='store_true', help='Check if Qpos interval matches Tpos interval when filtering blat results')
        parser.add_argument('--print-filtering-error', default=False, action='store_true', help='Display problematic data for blat filtering results')
        parser.add_argument('--time', default=False, action='store_true', help='Display elapsed time')
        parser.add_argument('--version', action='version', version=f'{prog}: {__version__}')
        self.args = parser.parse_args()

        self.__config(parser)

    def __config(self, parser):
        cfg_file = 'fuseq.cfg'
        cfg1 = f'os.path.abspath(".")/{cfg_file}'
        cfg2 = os.environ["FUSEQ_CFG"] if 'FUSEQ_CFG' in os.environ else ''
        cfg3 = os.path.expanduser(f'~/.config/{cfg_file}')

        # Check config file
        if os.path.exists(cfg1):
            cfg_path = cfg1
        elif os.path.exists(cfg2):
            cfg_path = cfg2
        elif os.path.exists(cfg3):
            cfg_path = cfg3
        else:
            return

        # Load config file
        cfg = configparser.ConfigParser()
        cfg.read(cfg_path)

        # Check section
        section = 'fuseq'
        if section not in cfg:
            return

        # Update self.args with a default value
        args_dic = vars(self.args)
        for key in cfg[section].keys():
            if key in args_dic:
                # Priority is an argument value rather than a config value
                if args_dic[key] != parser.get_default(key):
                    continue
                try:
                    if type(args_dic[key]) == bool:
                        value = cfg.getboolean(section, key)
                    elif isinstance(args_dic[key], int):
                        value = cfg.getint(section, key)
                    elif isinstance(args_dic[key], str):
                        value = cfg[section][key]
                    else:
                        print('[Error] at config function')
                        exit(1)
                    args_dic[key] = value
                except Exception:
                    print(f'[Error] Invalid value "{cfg[section][key]}" for "{key}" in config file')
                    exit(1)
            else:
                print(f'[Error] Invalid key "{key}" in config file')
                exit(1)

    def __check(self):
        args = self.args
        Checker.isdir(args.genomon_root_directory)
        Checker.isfile(args.reference)
        if args.shirokane:
            Checker.onshirokane()

    def __get_lines(self, line_str):
        """line_str='1,5,3,8-10,12' => [1, 3, 5, 8, 9, 10, 12]"""
        if not line_str:
            return []
        line_sp = line_str.split(',')
        lines = []
        for line in line_sp:
            m = re.match(r'^(\d+)-(\d+)$', line)
            if m:
                lines += list(range(int(m.group(1)), int(m.group(2)) + 1))
            elif not re.match(r'^\d+$', line):
                print(f'Not a line number: {line}')
                exit(1)
            else:
                lines.append(int(line))
        lines = sorted(list(set(lines)))
        # Add header line
        if lines[0] != 1:
            lines.insert(0, 1)
        return lines

    def __create(self):
        args = self.args
        num_coll_parallels, num_blat_parallels = self.__get_parallels(args)

        # Change required arguments
        args.fuseq_root_dir = os.path.abspath(args.fuseq_root_directory)
        args.genomon_root_dir = os.path.abspath(args.genomon_root_directory)
        # Change options
        args.blat_opts = self.__modify_blat_opts(args.blat_options)
        args.star_dir = os.path.abspath(args.star_directory) if args.star_directory else ''
        args.mf_lines = self.__get_lines(args.lines)
        args.fusion_file = os.path.abspath(args.fusion_file) if args.fusion_file else ''
        args.delete_work = not args.no_delete_work
        # restart_blat
        # restart_filter
        # collection_processes
        args.use_filt = False if args.no_use_filt else True
        args.readname_filt = args.readname
        args.seq_filt = args.sequence
        args.bp_start_extn = args.start
        args.bp_end_extn = args.end
        args.reference = os.path.abspath(args.reference)
        args.on_shirokane = args.shirokane
        args.num_coll_parallels = num_coll_parallels
        args.num_blat_parallels = num_blat_parallels
        args.check_pos_intvl = False if args.no_check_position_interval else True
        args.print_filt_err = args.print_filtering_error
        args.print_time = args.time

        # Delete
        del args.fuseq_root_directory, args.genomon_root_directory, \
        args.blat_options, args.star_directory, args.lines, args.no_delete_work, \
        args.collection_processes, args.no_use_filt, args.readname, args.sequence, \
        args.start, args.end, args.shirokane, args.collection_tasks, args.blat_tasks, \
        args.no_check_position_interval, args.print_filtering_error, args.time \

        # Add
        args.is_restart = True if args.restart_blat or args.restart_filter else False
        args.fuseq_filename = 'fusion_sequence.txt'
        args.work_dirname = 'work_restart'
        args.swork_dirname = 'splits'

    def __get_parallels(self, args):
        if args.shirokane:
            return args.collection_tasks, args.blat_tasks
        else:
            return args.collection_processes, 1

    def __modify_blat_opts(self, blat_opts):
        if not blat_opts:
            return ''
        blat_opts = [opt if opt[0] == '-' else f'-{opt}'
                     for opt in re.split(' +', blat_opts.strip(' '))]
        return ' '.join(blat_opts)

    def refer(self):
        return self.args

    def copy(self):
        return argparse.Namespace(**vars(self.args))  # Deep copy
