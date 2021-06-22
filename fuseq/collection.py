import csv
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
                breakinfo.append({
                    'linenr': linenr, 'sample': sample,
                    'chr1': chr1, 'bp1': bp1, 'strand1': strand1, 'gene1': gene1, 'junc1': junc1,
                    'chr2': chr2, 'bp2': bp2, 'strand2': strand2, 'gene2': gene2, 'junc2': junc2})
        return breakinfo

    def __create_script(self, breakinfo):
        # Commands for filtering
        readname_filt_cmd = \
            f'[ "$readname" != \'{self.params.readname_filt}\' ] && continue' \
            if self.params.readname_filt else ''
        seq_filt_cmd = \
            f'[ "$seq" != \'{self.params.seq_filt}\' ] && continue' \
            if self.params.seq_filt else ''

        # Commands
        cmd_head = '#!/bin/bash\n\nset -eu\n\n'
        cmd_main = '''\
chr1='{chr1}'
chr2='{chr2}'
bp1='{bp1}'
bp2='{bp2}'
jun_path='{jun_path}'
sam_path="${{jun_path%\\.*}}.sam"
out_path='{out_path}'

touch "$out_path"
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
\n
'''

        line_cnt = len(breakinfo)
        n_parallels = min(line_cnt, self.params.num_coll_parallels)
        width = len(str(n_parallels))
        # Determine the number of lines each process
        lines_each_proc = line_cnt // n_parallels
        n_plus1 = line_cnt - lines_each_proc * n_parallels
        if n_plus1 == 0:
            heads = [i * lines_each_proc for i in range(n_parallels + 1)]
        else:
            plus1lines_each_proc = lines_each_proc + 1
            total_plus1lines = plus1lines_each_proc * n_plus1
            n_plus0 = n_parallels - n_plus1
            heads = \
                [i * plus1lines_each_proc for i in range(n_plus1)] + \
                [total_plus1lines + i * lines_each_proc for i in range(n_plus0 + 1)]

        jun_dic = {}
        out_paths = []
        for i, (head, tail) in enumerate(zip(heads, heads[1:])):
            out_path = f'{self.params.swork_dir}/{self.out_file}{str(i+1).zfill(width)}'
            script_path = f'{out_path}.sh'
            out_paths.append(out_path)
            with open(script_path, 'w') as f:
                f.write(cmd_head.format(out_path=out_path))
                for d in breakinfo[head:tail]:
                    [linenr, sample, chr1, bp1, strand1, _, _,
                     chr2, bp2, strand2, _, _] = d.values()
                    if sample not in jun_dic:
                        jun_dic[sample] = glob.glob(f'{self.star_dir}/{sample}/*.junction')[0]
                    jun_path = jun_dic[sample]
                    bp1_arng = str(int(bp1) + 1) if strand1 == '+' else str(int(bp1) - 1)
                    bp2_arng = str(int(bp2) + 1) if strand2 == '+' else str(int(bp2) - 1)
                    cmd = cmd_main.format(linenr=linenr, chr1=chr1, bp1=bp1_arng, chr2=chr2, bp2=bp2_arng,
                                          jun_path=jun_path, out_path=out_path,
                                          readname_filt_cmd=readname_filt_cmd, seq_filt_cmd=seq_filt_cmd)
                    f.write(cmd)
            os.chmod(script_path, 0o0755)

        return out_paths

    def __task(self, i_proc, out_path, errs, rcs):
        cmd = f'bash {out_path}.sh'
        _, err, rc = self._run_cmd(cmd, 'collection', ignore_err=True)
        errs[i_proc] = err
        rcs[i_proc] = rc

    def __collect(self, breakinfo):
        """Collect data for Blat input"""
        # Create scripts
        out_paths = self.__create_script(breakinfo)
        n_parallels = len(out_paths)

        # Run scripts
        if self.params.on_shirokane:
            width = len(str(n_parallels))
            coll_path = out_paths[0][:-width]
            cmd = '''\
#!/usr/local/bin/nosh
#$ -S /usr/local/bin/nosh
#$ -cwd
#$ -l s_vmem=4G,mem_req=4G
#$ -e {coll_path}.log
#$ -o {coll_path}.log
id=$(printf "%0{width}d" ${{SGE_TASK_ID}})
bash {coll_path}${{id}}.sh
'''.format(coll_path=coll_path, width=width)
            script_path = f'{coll_path}.sh'
            self._run_cmd_on_uge(cmd, script_path, n_parallels, 'collection_uge')

        else:
            manager = multiprocessing.Manager()
            errs = manager.dict()
            rcs = manager.dict()
            jobs = []
            for i in range(n_parallels):
                p = multiprocessing.Process(
                    target=self.__task,
                    args=(i, out_paths[i], errs, rcs))
                jobs.append(p)
                p.start()
            for job in jobs:
                job.join()

            # Check return codes
            has_err = False
            for i in range(n_parallels):
                if rcs[i] != 0:
                    print('[Error] Return code is not 0 at collection script')
                    print(f'err: {errs[i]}, rc: {rcs[i]}')
                    has_err = True
            if has_err:
                exit(1)

        return out_paths

    def __add_count_to(self, breakinfo):
        coll_path = f'{self.params.work_dir}/{self.out_file}'
        cnts = [0] * len(breakinfo)
        with open(coll_path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            prev_cnt = 0
            tgt_linenr = 2  # first line is header line
            for row in reader:
                sp = row[0].split('_')[0].split('-')  # row[0]=2-1_READNAME => sp=[2,1]
                cur_linenr = int(sp[0][1:])           # sp[0][0] = '>'
                cur_cnt = int(sp[1])
                if cur_linenr != tgt_linenr:
                    cnts[tgt_linenr - 2] = prev_cnt
                    tgt_linenr = cur_linenr
                prev_cnt = cur_cnt
                next(reader)  # sequence data
            cnts[tgt_linenr - 2] = prev_cnt
        for i, cnt in enumerate(cnts):
            breakinfo[i]['cnt'] = cnt

    def __concat(self, inp_paths):
        inp_files = ' '.join([os.path.basename(path) for path in inp_paths])
        cmd = '''\
#!/bin/bash
set -eu
cd {swork_dir}
cat {inp_files} > ../{out_file}
'''.format(swork_dir=self.params.swork_dir, inp_files=inp_files, out_file=self.out_file)
        self._run_cmd(cmd, 'cat_coll_files')

    @Timer('collection')
    def run(self):
        self.__create_symlinks()
        breakinfo = self.__get_breakinfo()
        coll_out_paths = self.__collect(breakinfo)
        self.__concat(coll_out_paths)
        self.__add_count_to(breakinfo)
        return breakinfo
