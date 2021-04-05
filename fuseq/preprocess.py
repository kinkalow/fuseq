import os
import subprocess
import multiprocessing

class Preprocess():
    """
    Create input data for blat
    """

    def __init__(self, mf_dir, mf_path, jun_dic, opts):
        self.mf_dir = mf_dir    # directory name containing merge_fusionfusion
        self.mf_path = mf_path  # merge_fusionfusion path
        self.jun_dic = jun_dic  # {directory name containing junction file, junction path}
        self.out_base_dir = opts.out_base_dir
        self.parallel_num = opts.preprocess_parallel_num
        self.out_dir = f'{self.out_base_dir}/{self.mf_dir}/_work'
        #self._proc()

    def get_data(self):
        chrs = list(map(lambda x: str(x), range(23))) + list('XY')
        data = []
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
                data.append([sample, chr1, bp1, chr2, bp2])
        return data

    def process(self, data):
        os.makedirs(self.out_dir, exist_ok=True)

        with open(self.mf_path, 'r') as f:
            line_cnt = sum([1 for _ in f])
        parallel_num = min(line_cnt, self.parallel_num)
        data_num = line_cnt // parallel_num if line_cnt % parallel_num == 0 else line_cnt // parallel_num + 1
        heads = [i * data_num for i in range(parallel_num)] + [line_cnt]

        def task(i_proc):
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

#[ -f "$out_path" ] && rm "$out_path"

for readname in $(cat "$jun_path" | awk '{{ \\
  if ( ($1 == "'$chr1'" && $2 == "'$bp1'" && $4 == "'$chr2'" && $5 == "'$bp2'") || \\
       ($1 == "'$chr2'" && $2 == "'$bp2'" && $4 == "'$chr1'" && $5 == "'$bp1'")    \\
     ) print $10 }}'); do
    out=$(grep "^$readname" "$sam_path" | awk '{{ if ($9 == 0) print $10 }}')
    echo ">$readname\\n$out" >> "$out_path"
done
'''
            outs = []
            head = heads[i_proc]
            tail = heads[i_proc + 1]
            digit = len(str(parallel_num).split('.')[0])
            out_path = f'{self.out_dir}/{str(i_proc).zfill(digit)}'
            if os.path.exists(out_path):
                os.remove(out_path)

            for [sample, chr1, bp1, chr2, bp2] in data[head:tail]:
                jun_path = self.jun_dic[sample]
                cmd = cmd_template.format(chr1=chr1, bp1=bp1, chr2=chr2, bp2=bp2, jun_path=jun_path, out_path=out_path)
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out, err = p.communicate()
                if p.returncode != 0:
                    print(f'[Error] return code: {p.returncode} for {sample} {chr1} {bp1} {chr2} {bp2}')
                    print(err)
                    continue
                outs.append([out])

        import time
        start = time.time()
        jobs = []
        for i in range(parallel_num):
            p = multiprocessing.Process(target=task, args=(i,))
            jobs.append(p)
            p.start()
        for p in jobs:
            p.join()
        print(f"{time.time() - start:.3f}[s]")

    def process_slow(self, data):
        # Extract readnames from junction file
        for [sample, chr1, bp1, chr2, bp2] in data:
            jun_dic = self.jun_dic[sample]
            readnames = []
            with open(jun_dic, 'r') as f_jun:
                match = False
                for row_jun in f_jun:
                    r = row_jun.rstrip('\n').split('\t')
                    if (r[0] == chr1 and r[1] == bp1 and r[3] == chr2 and r[4] == bp2) or \
                       (r[3] == chr1 and r[4] == bp1 and r[0] == chr2 and r[1] == bp2):
                        readnames.append(r[9])
                        match = True
                if not match:
                    print(f'[Error] No data for {chr1} {bp1} {chr2} {bp2} in {jun_dic}')
                    exit(1)

            # Create data for blat input
            path_sam = jun_dic[:-8] + 'sam'
            results = []
            with open(path_sam, 'r') as f_sam:
                # Skip header lines
                while True:
                    pos = f_sam.tell()
                    row = f_sam.readline()
                    if row[0] != '@':
                        break
                # Alignment lines
                for readname in readnames:
                    f_sam.seek(pos)
                    match = False
                    for row_sam in f_sam:
                        if row_sam.startswith(readname):
                            r = row_sam.rstrip('\n').split('\t')
                            if r[8] == '0':
                                results.append([f'>{readname}\n{r[9]}'])
                                match = True
                    if not match:
                        print(f'[Error] No data for {readname} in {path_sam}')
                        exit(1)
                print(results)
                exit()
