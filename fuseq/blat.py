import json
import os
import subprocess
import multiprocessing
import shutil
from fuseq.timer import Timer
from fuseq.checker import Checker

class Blat():

    def __init__(self, mf_path, jun_dic, work_dir, fuseq_path, params):
        self.mf_path = mf_path  # merge_fusionfusion path
        self.jun_dic = jun_dic  # jun_dic = {directory name containing junction file, full path of junction file}
        self.work_dir = work_dir
        self.fuseq_path = fuseq_path
        # Parameters
        self.params = params
        # File names / Full path
        self.files = {'params': 'params', 'coll': 'collect',
                      'breakinfo': 'breakinfo', 'blat': 'blat',
                      'filtsome': 'filter_some', 'filtemp': 'filter_empty'}
        self.breakinfo_path = f'{self.work_dir}/{self.files["breakinfo"]}'
        # Timer
        self.print_time = True
        self.time_coll = Timer('Collection')
        self.time_blat = Timer('Blat')
        self.time_filter = Timer('Filter')

    #
    # Utility
    #

    def save_params(self):
        path = f'{self.work_dir}/{self.files["params"]}'
        with open(path, 'w') as f:
            d = vars(self.params)
            max_len = max([len(k) for k in d.keys()])
            for key, val in sorted(d.items(), key=lambda x: x[0]):
                space = ' ' * (max_len - len(key))
                f.write(f'{key}{space}: {val}\n')

    def save_breakinfo(self, breakinfo):
        with open(self.breakinfo_path, 'w') as f:
            json.dump(breakinfo, f)

    def load_breakinfo(self):
        Checker.isfile(self.breakinfo_path)
        with open(self.breakinfo_path) as f:
            breakinfo = json.load(f)
        return breakinfo

    def delete_work_dir(self, make_work_dir=False):
        shutil.rmtree(self.work_dir, ignore_errors=True)
        if make_work_dir:
            os.makedirs(self.work_dir)

    #
    # Main processing functions
    #

    def get_breakinfo(self):
        chrs = list(map(lambda x: str(x), range(23))) + list('XY')
        breakinfo = []
        with open(self.mf_path, 'r') as f_mf:
            # Skip header line
            f_mf.readline()
            # Main lines
            for row_mf in f_mf:
                # Extract chrs and bps from merge_fusionfusion file
                row_mf = row_mf.rstrip('\n').split('\t')
                if row_mf[1] not in chrs or row_mf[4] not in chrs:
                    continue
                [sample, chr1, bp1, strand1, chr2, bp2, strand2] = row_mf[0:7]
                bp1 = str(int(bp1) + 1) if strand1 == '+' else str(int(bp1) - 1)
                bp2 = str(int(bp2) + 1) if strand2 == '+' else str(int(bp2) - 1)
                breakinfo.append([sample, chr1, bp1, strand1, chr2, bp2, strand2])
        return breakinfo

    def collect(self, breakinfo):
        """Collect data for Blat input"""
        self.time_coll.start()

        with open(self.mf_path, 'r') as f:
            line_cnt = sum([1 for _ in f])
        coll_procs = min(line_cnt, self.params.coll_procs)
        data_num = line_cnt // coll_procs if line_cnt % coll_procs == 0 else line_cnt // coll_procs + 1
        heads = [i * data_num for i in range(coll_procs)] + [line_cnt]

        def task(i_proc, sender):
            cmd_template = '''\
#!/bin/bash

set -eu

chr1='{chr1}'
chr2='{chr2}'
bp1='{bp1}'
bp2='{bp2}'
jun_path='{jun_path}'
sam_path="${{jun_path%\\.*}}.sam"
out_path='{out_path}'

cnt='0'
for readname in $(cat "$jun_path" | awk '{{ \\
  if ( ($1 == "'$chr1'" && $2 == "'$bp1'" && $4 == "'$chr2'" && $5 == "'$bp2'") || \\
       ($1 == "'$chr2'" && $2 == "'$bp2'" && $4 == "'$chr1'" && $5 == "'$bp1'")    \\
     ) print $10 }}'); do
    seqs=$(grep "^$readname" "$sam_path" | awk '{{ if ($7 != "=" && $9 == 0 && $15 != "XS:A:+") print $10 }}')
    [ -z "$seqs" ] && continue
    for seq in $seqs; do
      cnt=$((cnt+1))
      printf ">{line}-${{cnt}}_$readname\\n$seq\\n" >> "$out_path"
    done
done
echo -n "$cnt"
'''
            results = []
            head = heads[i_proc]
            tail = heads[i_proc + 1]
            digit = len(str(coll_procs).split('.')[0])
            out_path = f'{self.work_dir}/{str(i_proc).zfill(digit)}'
            if os.path.exists(out_path):
                os.remove(out_path)

            for step, [sample, chr1, bp1, strand1, chr2, bp2, strand2] in enumerate(breakinfo[head:tail]):
                line = head + step + 1
                jun_path = self.jun_dic[sample]
                cmd = cmd_template.format(line=line, chr1=chr1, bp1=bp1, chr2=chr2, bp2=bp2,
                                          jun_path=jun_path, out_path=out_path)
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out, err = p.communicate()
                results.append({'line': line, 'ret': p.returncode,
                                'err': err.decode(), 'cnt': out.decode(),
                                'sample': sample,
                                'chr1': chr1, 'bp1': bp1, 'strand1': strand1,
                                'chr2': chr2, 'bp2': bp2, 'strand2': strand2})
            sender.send(results)
            sender.close()

        jobs = []
        pipes = []
        for i in range(coll_procs):
            receiver, sender = multiprocessing.Pipe(False)
            p = multiprocessing.Process(target=task, args=(i, sender))
            jobs.append(p)
            pipes.append(receiver)
            p.start()

        # results
        breakinfo = []
        has_err = False
        for pipe in pipes:
            results = pipe.recv()  # Wait for the results to be sent
            for r in results:
                if r['ret'] != 0 or len(r['err']) != 0 or not r['cnt'].isdecimal():
                    print(f'[Error] {str(r)}')
                    has_err = True
                r.pop('ret')
                r.pop('err')
                breakinfo.append(r)

        for job in jobs:
            job.join()

        if has_err:
            exit(1)

        self.time_coll.end()
        if self.print_time:
            self.time_coll.print()

        return breakinfo

    def concat(self):
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
cat {inp_files} > {out_file}
'''.format(work_dir=self.work_dir,
           inp_files=' '.join(map(str, range(self.params.coll_procs))),
           out_file=os.path.basename(self.files['coll']))
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
        _, err = p.communicate()
        if p.returncode != 0:
            print('[Error] at cat runtime')
            print(err)
            exit(1)

    def blat(self):
        self.time_blat.start()

        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
blat {blat_opt} -noHead {reference} {inp_file} {out_file}
#blat -noHead {reference} {inp_file} {out_file}
'''.format(work_dir=self.work_dir, blat_opt=self.params.blat_opt, reference=self.params.reference,
           inp_file=self.files['coll'], out_file=self.files['blat'])
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
        _, err = p.communicate()
        if p.returncode != 0:
            print('[Error] at blat runtime')
            print(err)
            exit(1)

        self.time_blat.end()
        if self.print_time:
            self.time_blat.print()

    # [1] bp1 and bp2 are in the range of Tstart and Tend
    # [2] chr1 and chr2 are the same as Tname
    # [3] No range if only one item satisfies both [1] and [2]
    # [4] In the case of three or more items that satisfy [1] and [2]
    #     when the range of [Qstart2,Qend2] is wider than that of [Qstart1,Qend1], [Qstart2,Qend2] is given priority
    #     [Qstart1,Qend1] is not displayed
    def blat_filter(self, breakinfo):
        self.time_filter.start()

        # Get read names
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
grep '^>' {inp_file} |  sed -e 's/^>//'
'''.format(work_dir=self.work_dir, inp_file=self.files['coll'])
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        readnames, err = p.communicate()
        if p.returncode != 0:
            print('[Error] at grep runtime')
            print(err)
            exit(1)
        readnames = readnames.decode().strip('\n').split('\n')

        # Check
        cnts = [int(d.get('cnt')) for d in breakinfo]
        assert(len(readnames) == sum(cnts))

        # Create a dictionary with readname in key and breakinfo index in value
        i_rn = 0
        rn2idx = dict()
        for i, cnt in enumerate(cnts):
            for _ in range(cnt):
                rn2idx[readnames[i_rn]] = i
                i_rn += 1

        def update_chr_bp(i, d):
            return d[i]['chr1'], d[i]['chr2'], int(d[i]['bp1']), int(d[i]['bp2'])

        def write_to_file(fws, breakinfo, readname, seq, poses_filt, other_info):
            chr1 = breakinfo['chr1']
            bp1 = breakinfo['bp1']
            strand1 = breakinfo['strand1']
            chr2 = breakinfo['chr2']
            bp2 = breakinfo['bp2']
            strand2 = breakinfo['strand2']
            mfline = breakinfo['line']
            readname = readname[readname.index('_') + 1:]  # Remove the leading unique number

            w1 = readname + '\n'
            w2 = ' '.join((chr1, bp1, strand1, chr2, bp2, strand2, f'mfline={mfline}')) + '\n'
            if poses_filt:
                fw = fws[0]
                w3 = ''
                for s, e, _ in poses_filt:
                    w3 += ' ' * (s - 1) + seq[s - 1:e] + ' ' + str([s, e]).replace(' ', '') + '\n'
            else:
                fw = fws[1]
                w3 = seq + '\n'
            w = w1 + w2 + w3
            #w += str(other_info) + '\n'
            fw.write(w + '\n')
            #print(w, end='')
            #print('-' * 126)

        # Filter base on qstat-qend range
        def pos_filter(poses):
            if len(poses) < 2:
                return []

            starts = [r[0] for r in poses]
            idx_asc = sorted(range(len(starts)), key=lambda k: starts[k])
            poses_asc = [poses[i] for i in idx_asc]  # poses[i] = (start, end, one_or_two)

            cur_idx = 0
            max_idx = len(idx_asc) - 1
            results = [0] * len(idx_asc)  # if 0/1, unnecessary/necessary elements
            while cur_idx < max_idx:
                s, e, _ = poses_asc[cur_idx]
                # Extract a location with the same start position and the largest end position
                for i in range(cur_idx + 1, max_idx + 1):
                    if poses_asc[i][0] != s:
                        i -= 1
                        break
                    if poses_asc[i][1] > e:
                        e = poses_asc[i][1]
                else:
                    results[cur_idx] = 1
                    break
                results[i] = 1
                # Find the next candidate
                for i in range(i + 1, max_idx + 1):
                    if poses_asc[i][1] > e:
                        break
                else:
                    break
                if i == max_idx:
                    results[i] = 1
                    break
                cur_idx = i

            poses_filt = [poses[i] for i, zero_or_one in enumerate(results) if zero_or_one == 1]
            if len(poses_filt) < 2:
                return []

            return poses_filt

        def pos_sort(poses_filt, breakinfo):
            if not poses_filt:
                return poses_filt, breakinfo

            # Check
            one_or_two_list = [i for _, _, i in poses_filt]
            if len(set(one_or_two_list)) != 2:
                print('[Error] poses_filt has some problem')
                print(poses_filt)
                exit(1)

            # Sort position and break information
            # if one_or_two == 1, sequence1 is displayed first before sequence2
            # if one_or_two == 2, sequence2 is displayed first before sequence1
            starts = [s for s, _, _ in poses_filt]
            idx_min = starts.index(min(starts))
            is_one_first = True if poses_filt[idx_min][2] == 1 else False
            if is_one_first:
                idxes = [i for i, one_or_two in enumerate(one_or_two_list) if one_or_two == 1] + \
                        [i for i, one_or_two in enumerate(one_or_two_list) if one_or_two == 2]
            else:
                idxes = [i for i, one_or_two in enumerate(one_or_two_list) if one_or_two == 2] + \
                        [i for i, one_or_two in enumerate(one_or_two_list) if one_or_two == 1]
            new_poses_filt = [poses_filt[idx] for idx in idxes]

            if poses_filt[0][2] == 1:
                new_breakinfo = breakinfo
            else:  # Reverse
                b = breakinfo.copy()
                b['chr2'], b['bp2'], b['strand2'], b['chr1'], b['bp1'], b['strand1'] = \
                b['chr1'], b['bp1'], b['strand1'], b['chr2'], b['bp2'], b['strand2']
                new_breakinfo = b

            return new_poses_filt, new_breakinfo

        def filter_and_write(fws, breakinfo, readname, seq, poses, other_info):
            poses_filt = pos_filter(poses)
            poses_filt, breakinfo = pos_sort(poses_filt, breakinfo)
            write_to_file(fws, breakinfo, readname, seq, poses_filt, other_info)
            poses.clear()

        # Full path
        blat_path = f'{self.work_dir}/{self.files["blat"]}'
        coll_path = f'{self.work_dir}/{self.files["coll"]}'
        filtsome_path = f'{self.work_dir}/{self.files["filtsome"]}'
        filtemp_path = f'{self.work_dir}/{self.files["filtemp"]}'
        fuseq_path = self.fuseq_path

        # Open
        fr_preproc = open(coll_path, 'r')
        fr_blat = open(blat_path, 'r')
        fw_filtsome = open(filtsome_path, 'w')
        fw_filtemp = open(filtemp_path, 'w')
        fw_filts = [fw_filtsome, fw_filtemp]

        # Filter and write
        is_current = True
        for row in fr_preproc:
            # Target read
            cur_readname = row.rstrip('\n')[1:]
            cur_seq = fr_preproc.readline().rstrip('\n')
            cur_chr1, cur_chr2, cur_bp1, cur_bp2 = update_chr_bp(rn2idx[cur_readname], breakinfo)
            cur_breakinfo = breakinfo[rn2idx[cur_readname]]

            poses = []
            other_info = []
            while True:
                if is_current:
                    s = fr_blat.readline().rstrip('\n').split('\t')
                    if s == ['']:  # Check if file pointer is last
                        break
                    readname = s[9]             # qname
                    pos_start = int(s[11]) + 1  # qstart
                    pos_end = int(s[12])        # qend
                    chr = s[13]                 # tname
                    bp_start = int(s[15])       # tstart
                    bp_end = int(s[16]) + 1     # tend ... NOTE: +1 increases the number of matches to Genomon results
                    assert(bp_start <= bp_end)
                if readname == cur_readname:
                    # Filter based on chr and tstart-tend range
                    if chr == cur_chr1 and bp_start <= cur_bp1 <= bp_end:
                        poses.append((pos_start, pos_end, 1))
                    elif chr == cur_chr2 and bp_start <= cur_bp2 <= bp_end:
                        poses.append((pos_start, pos_end, 2))
                    if chr == cur_chr1:
                        other_info.append((bp_start, bp_end, 1))
                    elif chr == cur_chr2:
                        other_info.append((bp_start, bp_end, 2))
                    is_current = True
                else:
                    # Write one read
                    filter_and_write(fw_filts, cur_breakinfo, cur_readname, cur_seq, poses, other_info)
                    is_current = False
                    break
        if is_current:
            # Write the remaining read
            filter_and_write(fw_filts, cur_breakinfo, cur_readname, cur_seq, poses, other_info)

        # Close
        fr_preproc.close()
        fr_blat.close()
        fw_filtsome.close()
        fw_filtemp.close()

        # Concat filtsome and filtemp
        filtsome_exists = os.path.exists(filtsome_path)
        filtemp_exists = os.path.exists(filtemp_path)
        if filtsome_exists and filtemp_exists:
            cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
cat {filtsome} {filtemp} > {fuseq}
'''.format(work_dir=self.work_dir, fuseq=fuseq_path, filtsome=filtsome_path, filtemp=filtemp_path)
            p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
            _, err = p.communicate()
            if p.returncode != 0:
                print('[Error] at filter runtime')
                print(err)
                exit(1)
        elif filtsome_exists:
            shutil.move(filtsome_path, fuseq_path)
        elif filtemp_exists:
            shutil.move(filtemp_path, fuseq_path)

        self.time_filter.end()
        if self.print_time:
            self.time_filter.print()

    #
    # Start functions
    #

    def run(self):
        # Clear and Save
        self.delete_work_dir(True)
        self.save_params()
        # Collect
        breakinfo = self.get_breakinfo()
        breakinfo = self.collect(breakinfo)
        self.concat()
        self.save_breakinfo(breakinfo)
        # Blat
        self.blat()
        # Filter
        self.blat_filter(breakinfo)
        # Delete
        if self.params.delete_work:
            self.delete_work_dir()

    def restart(self):
        breakinfo = self.load_breakinfo()
        self.save_params()
        if self.params.restart_blat:
            self.blat()
        self.blat_filter(breakinfo)
        if self.params.delete_work:
            self.delete_work_dir()
