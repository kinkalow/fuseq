import os
import subprocess
import multiprocessing
import shutil
from fuseq.timer import Timer

class Blat():

    def __init__(self, mf_path, jun_dic, work_dir, fuseq_path, opts):
        #self.mf_dir = mf_dir    # directory name containing merge_fusionfusion
        self.mf_path = mf_path  # merge_fusionfusion path
        self.jun_dic = jun_dic  # {directory name containing junction file, junction path}
        self.parallel_num = opts.preprocess_parallel_num
        self.work_dir = work_dir
        self.reference = opts.reference
        self.files = {'preproc': 'preprocess', 'blat': 'blat',
                      'filtsome': 'filter_some', 'filtemp': 'filter_empty',
                      'fuseq': fuseq_path}
        self.time_preproc = Timer('Preprocess')
        self.time_blat = Timer('Blat')
        self.time_filter = Timer('Filter')

    def get_data(self):
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

    def clear_work_dir(self, make_work_dir=False):
        shutil.rmtree(self.work_dir, ignore_errors=True)
        if make_work_dir:
            os.makedirs(self.work_dir)

    def preprocess(self, breakinfo):
        with open(self.mf_path, 'r') as f:
            line_cnt = sum([1 for _ in f])
        parallel_num = min(line_cnt, self.parallel_num)
        data_num = line_cnt // parallel_num if line_cnt % parallel_num == 0 else line_cnt // parallel_num + 1
        heads = [i * data_num for i in range(parallel_num)] + [line_cnt]

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
    out=$(grep "^$readname" "$sam_path" | awk '{{ if ($7 != "=" && $9 == 0 && $15 != "XS:A:+") print $10 }}')
    [ -z "$out" ] && continue
    cnt=$((cnt+1))
    echo ">$readname\\n$out" >> "$out_path"
done
echo -n "$cnt"
'''
            results = []
            head = heads[i_proc]
            tail = heads[i_proc + 1]
            digit = len(str(parallel_num).split('.')[0])
            out_path = f'{self.work_dir}/{str(i_proc).zfill(digit)}'
            if os.path.exists(out_path):
                os.remove(out_path)

            for step, [sample, chr1, bp1, strand1, chr2, bp2, strand2] in enumerate(breakinfo[head:tail]):
                jun_path = self.jun_dic[sample]
                cmd = cmd_template.format(chr1=chr1, bp1=bp1, chr2=chr2, bp2=bp2,
                                          jun_path=jun_path, out_path=out_path)
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out, err = p.communicate()
                results.append({'line': head + step + 1, 'ret': p.returncode,
                                'err': err.decode(), 'cnt': out.decode(),
                                'sample': sample,
                                'chr1': chr1, 'bp1': bp1, 'strand1': strand1,
                                'chr2': chr2, 'bp2': bp2, 'strand2': strand2})
            sender.send(results)
            sender.close()

        jobs = []
        pipes = []
        for i in range(parallel_num):
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

        return breakinfo

    def cat(self):
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
cat {inp_files} > {out_file}
'''.format(work_dir=self.work_dir,
           inp_files=' '.join(map(str, range(self.parallel_num))),
           out_file=self.files['preproc'])
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
        p.communicate()
        if p.returncode != 0:
            print('[Error] at cat runtime')

    def blat(self):
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
blat -noHead {reference} {inp_file} {out_file}
'''.format(work_dir=self.work_dir, reference=self.reference,
           inp_file=self.files['preproc'], out_file=self.files['blat'])
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        if p.returncode != 0:
            print('[Error] at blat runtime')

    # [1] bp1 and bp2 are in the range of Tstart and Tend
    # [2] chr1 and chr2 are the same as Tname
    # [3] No range if only one item satisfies both [1] and [2]
    # [4] In the case of three or more items that satisfy [1] and [2]
    #     when the range of [Qstart2,Qend2] is wider than that of [Qstart1,Qend1], [Qstart2,Qend2] is given priority
    #     [Qstart1,Qend1] is not displayed
    def blat_filter(self, breakinfo):
        # Get read names
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
grep '^>' {inp_file} |  sed -e 's/^>//'
'''.format(work_dir=self.work_dir, inp_file=self.files['preproc'])
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        readnames, err = p.communicate()
        if p.returncode != 0:
            print('[Error] at grep runtime')
            print(err)
            exit(1)
        readnames = readnames.decode().strip('\n').split('\n')

        cnts = [int(d.get('cnt')) for d in breakinfo]
        assert(len(readnames) == sum(cnts))

        # Create a dictionary that points to the index of data containing chr and bp from the read name
        i_rn = 0
        rn2idx = dict()
        for i, cnt in enumerate(cnts):
            for _ in range(cnt):
                rn2idx[readnames[i_rn]] = i
                i_rn += 1

        def update_chr_bp(i, d):
            return d[i]['chr1'], d[i]['chr2'], int(d[i]['bp1']), int(d[i]['bp2'])

        def write_to_file(fws, breakinfo, readname, seq, poses_filt):
            chr1 = breakinfo['chr1']
            bp1 = breakinfo['bp1']
            strand1 = breakinfo['strand1']
            chr2 = breakinfo['chr2']
            bp2 = breakinfo['bp2']
            strand2 = breakinfo['strand2']
            mfline = breakinfo['line']
            w1 = readname + '\n'
            w2 = ' '.join((chr1, bp1, strand1, chr2, bp2, strand2, f'mfline={mfline}')) + '\n'
            if poses_filt:
                fw = fws[0]
                w3 = ''
                for s, e in poses_filt:
                    w3 += ' ' * (s - 1) + seq[s - 1:e] + ' ' + str([s, e]).replace(' ', '') + '\n'
            else:
                fw = fws[1]
                w3 = seq + '\n'
            w = w1 + w2 + w3
            fw.write(w + '\n')
            #print(w, end='')
            #print('-' * 126)

        # Filter base on qstat-qend range
        def pos_filter(poses):
            if len(poses) < 2:
                return []

            starts = [r[0] for r in poses]
            idx_asc = sorted(range(len(starts)), key=lambda k: starts[k])
            poses_asc = [poses[i] for i in idx_asc]  # poses[i] = (start, end)

            cur_idx = 0
            max_idx = len(idx_asc) - 1
            results = [0] * len(idx_asc)  # if 0/1, unnecessary/necessary elements
            while cur_idx < max_idx:
                s, e = poses_asc[cur_idx]
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

        def filter_and_write(fws, breakinfo, readname, seq, poses):
            poses_filt = pos_filter(poses)
            write_to_file(fws, breakinfo, readname, seq, poses_filt)
            poses.clear()

        # File information
        file_blat = f'{self.work_dir}/{self.files["blat"]}'
        file_preproc = f'{self.work_dir}/{self.files["preproc"]}'
        file_filtsome = f'{self.work_dir}/{self.files["filtsome"]}'
        file_filtemp = f'{self.work_dir}/{self.files["filtemp"]}'
        file_fuseq = self.files['fuseq']

        # Open
        fr_preproc = open(file_preproc, 'r')
        fr_blat = open(file_blat, 'r')
        fw_filtsome = open(file_filtsome, 'w')
        fw_filtemp = open(file_filtemp, 'w')
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
            while True:
                if is_current:
                    s = fr_blat.readline().rstrip('\n').split('\t')
                    if s == ['']:
                        break
                    readname = s[9]             # qname
                    pos_start = int(s[11]) + 1  # qstart
                    pos_end = int(s[12])        # qend
                    chr = s[13]                 # tname
                    bp_start = int(s[15])       # tstart
                    bp_end = int(s[16])         # tend
                if readname == cur_readname:
                    # Filter based on chr and tstart-tend range
                    if (chr == cur_chr1 and bp_start <= cur_bp1 <= bp_end) or \
                       (chr == cur_chr2 and bp_start <= cur_bp2 <= bp_end):
                        poses.append((pos_start, pos_end))
                    is_current = True
                else:
                    # Write one read
                    filter_and_write(fw_filts, cur_breakinfo, cur_readname, cur_seq, poses)
                    is_current = False
                    break
        if is_current:
            # Write the remaining read
            filter_and_write(fw_filts, cur_breakinfo, cur_readname, cur_seq, poses)

        # Close
        fr_preproc.close()
        fr_blat.close()
        fw_filtsome.close()
        fw_filtemp.close()

        # Concat filtsome and filtemp
        filtsome_exists = os.path.exists(file_filtsome)
        filtemp_exists = os.path.exists(file_filtemp)
        if filtsome_exists and filtemp_exists:
            cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
cat {filtsome} {filtemp} > {fuseq}
'''.format(work_dir=self.work_dir, fuseq=file_fuseq, filtsome=file_filtsome, filtemp=file_filtemp)
            p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
            _, err = p.communicate()
            if p.returncode != 0:
                print('[Error] at filter runtime')
                print(err)
                exit(1)
        elif filtsome_exists:
            shutil.move(file_filtsome, file_fuseq)
        elif filtemp_exists:
            shutil.move(file_filtemp, file_fuseq)

    def run(self):
        breakinfo = self.get_data()
        self.clear_work_dir(True)

        self.time_preproc.start()
        breakinfo = self.preprocess(breakinfo)
        self.cat()
        self.time_preproc.print()

        self.time_blat.start()
        self.blat()
        self.time_blat.print()

        #import json
        #with open(f'{self.work_dir}/breakinfo', 'w') as f:
        #    json.dump(breakinfo, f)
        #import json
        #with open(f'{self.work_dir}/breakinfo') as f:
        #    breakinfo = json.load(f)
        #breakinfo = [{'line': 1, 'cnt': '4', 'sample': 'RNA_PVS-CC001-01_20160520_01', 'chr1': '2', 'bp1': '216261859', 'strand1': '-', 'chr2': '9', 'bp2': '131833825', 'strand2': '+'}, {'line': 2, 'cnt': '20', 'sample': 'RNA_PVS-CC001-01_20160520_01', 'chr1': '3', 'bp1': '197592293', 'strand1': '-', 'chr2': '3', 'bp2': '197602647', 'strand2': '+'}, {'line': 3, 'cnt': '69', 'sample': 'RNA_PVS-CC001-01_20160520_01', 'chr1': '1', 'bp1': '110464617', 'strand1': '+', 'chr2': '2', 'bp2': '216259443', 'strand2': '+'}, {'line': 4, 'cnt': '6', 'sample': 'RNA_PVS-CC001-01_20160520_01', 'chr1': '1', 'bp1': '219366423', 'strand1': '-', 'chr2': '1', 'bp2': '219414651', 'strand2': '+'}, {'line': 5, 'cnt': '4', 'sample': 'RNA_PVS-CC001-01_20160520_01', 'chr1': '12', 'bp1': '29492783', 'strand1': '-', 'chr2': '12', 'bp2': '29494152', 'strand2': '+'}, {'line': 6, 'cnt': '6', 'sample': 'RNA_PVS-CC001-01_20160520_01', 'chr1': '14', 'bp1': '22111807', 'strand1': '+', 'chr2': '14', 'bp2': '22977590', 'strand2': '-'}, {'line': 7, 'cnt': '5', 'sample': 'RNA_PVS-CC001-01_20160520_01', 'chr1': '12', 'bp1': '29492783', 'strand1': '-', 'chr2': '12', 'bp2': '29494750', 'strand2': '+'}]

        self.time_filter.start()
        self.blat_filter(breakinfo)
        self.time_filter.print()

        self.clear_work_dir()
