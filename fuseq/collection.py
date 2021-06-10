import glob
import multiprocessing
import os
from fuseq.timer import Timer
from fuseq.base import Base


class Collection(Base):

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.input_dir = f'{os.path.dirname(params.work_dir)}/input'
        self.mf_path = f'{self.input_dir}/fusion.txt'
        self.star_dir = f'{self.input_dir}/{os.path.basename(params.inputs["star_dir"])}'
        self.out_file = self.files['coll']

    def __create_symlinks(self):
        # Create input directory
        os.makedirs(self.input_dir, exist_ok=True)
        # Remove files
        if os.path.lexists(self.mf_path):
            os.remove(self.mf_path)
        if os.path.lexists(self.star_dir):
            os.remove(self.star_dir)
        # Create star symbolic file
        os.symlink(self.params.inputs['star_dir'], self.star_dir)
        # Create fusion symbolic file or new file
        if not self.params.mf_lines:
            os.symlink(self.params.inputs['mf_path'], self.mf_path)
        else:
            # Extract specified lines from a fusion file and write them
            idx = 0
            with open(self.params.inputs['mf_path'], 'r') as fr:
                with open(self.mf_path, 'w') as fw:
                    for i, row in enumerate(fr, start=1):
                        if i == self.params.mf_lines[idx]:
                            fw.write(row)
                            idx += 1
                            if idx > len(self.params.mf_lines) - 1:
                                break

    def __get_breakinfo(self):
        chrs = list(map(lambda x: str(x), range(23))) + list('XY')
        breakinfo = []
        with open(self.mf_path, 'r') as f_mf:
            # Skip header line
            f_mf.readline()
            # Main lines
            linenr = 1
            for row_mf in f_mf:
                # Extract chrs and bps from merge_fusionfusion file
                linenr += 1
                row_mf = row_mf.rstrip('\n').split('\t')
                if row_mf[1] not in chrs or row_mf[4] not in chrs:
                    continue
                [sample, chr1, bp1, strand1, chr2, bp2, strand2] = row_mf[0:7]
                [gene1, junc1, gene2, junc2] = row_mf[8:12]
                bp1 = str(int(bp1) + 1) if strand1 == '+' else str(int(bp1) - 1)
                bp2 = str(int(bp2) + 1) if strand2 == '+' else str(int(bp2) - 1)
                breakinfo.append([linenr, sample, chr1, bp1, strand1, gene1, junc1, chr2, bp2, strand2, gene2, junc2])
        return breakinfo

    def __collect(self, breakinfo):
        """Collect data for Blat input"""
        with open(self.mf_path, 'r') as f:
            line_cnt = sum([1 for _ in f])
        coll_procs = min(line_cnt, self.params.collection_process)
        data_num = line_cnt // coll_procs if line_cnt % coll_procs == 0 else line_cnt // coll_procs + 1
        heads = [i * data_num for i in range(coll_procs)] + [line_cnt]

        # Filter commands
        readname_filt_cmd = ''
        if self.params.readname_filt:
            readname_filt_cmd = f'[ "$readname" != \'{self.params.readname_filt}\' ] && continue'
        seq_filt_cmd = ''
        if self.params.seq_filt:
            seq_filt_cmd = f'[ "$seq" != \'{self.params.seq_filt}\' ] && continue'

        jobs = []
        pipes = []
        for i in range(coll_procs):
            receiver, sender = multiprocessing.Pipe(False)
            p = multiprocessing.Process(
                target=self.__task,
                args=(i, sender, breakinfo, coll_procs, heads,
                      readname_filt_cmd, seq_filt_cmd))
            jobs.append(p)
            pipes.append(receiver)
            p.start()

        # results
        out_paths = []
        breakinfo = []
        has_err = False
        for pipe in pipes:
            out_path, results = pipe.recv()  # Wait for the results to be sent
            for r in results:
                if r['ret'] != 0 or len(r['err']) != 0 or not r['cnt'].isdecimal():
                    print(f'[Error] {str(r)}')
                    has_err = True
                r.pop('ret')
                r.pop('err')
                breakinfo.append(r)
            if os.path.exists(out_path):
                out_paths.append(out_path)

        for job in jobs:
            job.join()

        if has_err:
            exit(1)

        if not out_paths:
            print('[Error] No input data for blat')
            exit(1)

        return out_paths, breakinfo

    def __task(self, i_proc, sender, breakinfo, coll_procs, heads, readname_filt_cmd, seq_filt_cmd):
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
    {readname_filt_cmd}
    seqs=$(grep "^$readname" "$sam_path" | awk '{{ if ($7 != "=" && $9 == 0 && $15 != "XS:A:+") print $10 }}')
    [ -z "$seqs" ] && continue
    for seq in $seqs; do
      {seq_filt_cmd}
      cnt=$((cnt+1))
      printf ">{linenr}-${{cnt}}_$readname\\n$seq\\n" >> "$out_path"
    done
done
echo -n "$cnt"
'''
        results = []
        head = heads[i_proc]
        tail = heads[i_proc + 1]
        digit = len(str(coll_procs).split('.')[0])
        out_path = f'{self.params.work_dir}/{self.out_file}{str(i_proc).zfill(digit)}'
        if os.path.exists(out_path):
            os.remove(out_path)

        jun_dic = {}
        for step, [linenr, sample, chr1, bp1, strand1, gene1, junc1,
                   chr2, bp2, strand2, gene2, junc2] in enumerate(breakinfo[head:tail]):
            if sample not in jun_dic:
                jun_dic[sample] = glob.glob(f'{self.star_dir}/{sample}/*.junction')[0]
            jun_path = jun_dic[sample]
            cmd = cmd_template.format(linenr=linenr, chr1=chr1, bp1=bp1, chr2=chr2, bp2=bp2,
                                      jun_path=jun_path, out_path=out_path,
                                      readname_filt_cmd=readname_filt_cmd, seq_filt_cmd=seq_filt_cmd)
            out, err, ret = self._run_cmd(cmd, ignore_err=True)
            results.append({'linenr': linenr, 'ret': ret, 'err': err, 'cnt': out, 'sample': sample,
                            'chr1': chr1, 'bp1': bp1, 'strand1': strand1, 'gene1': gene1, 'junc1': junc1,
                            'chr2': chr2, 'bp2': bp2, 'strand2': strand2, 'gene2': gene2, 'junc2': junc2})
        sender.send((out_path, results))
        sender.close()

    def __concat(self, inp_paths):
        inp_files = ' '.join([os.path.basename(path) for path in inp_paths])
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
cat {inp_files} > {out_file}
'''.format(work_dir=self.params.work_dir, inp_files=inp_files, out_file=self.out_file)
        self._run_cmd(cmd, 'cat_coll_files')

    @Timer('collection')
    def run(self):
        self.__create_symlinks()
        breakinfo = self.__get_breakinfo()
        coll_out_paths, breakinfo = self.__collect(breakinfo)
        self.__concat(coll_out_paths)
        return breakinfo
