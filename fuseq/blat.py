import os
import subprocess
import glob
import time
import multiprocessing

class Blat():
    """
    Create input data for blat
    """

    def __init__(self, mf_path, jun_dic, work_dir, opts):
        #self.mf_dir = mf_dir    # directory name containing merge_fusionfusion
        self.mf_path = mf_path  # merge_fusionfusion path
        self.jun_dic = jun_dic  # {directory name containing junction file, junction path}
        self.parallel_num = opts.preprocess_parallel_num
        self.work_dir = work_dir
        self.reference = opts.reference
        self.elapsed_time = {'preprocess': 0, 'blat': 0, 'filter': 0}
        self.files = {'preproc': 'preprocess', 'blat': 'blat.org', 'filter': 'blat'}

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

    def preprocess(self, data):
        for file in glob.glob(self.work_dir + '/*'):
            os.remove(file)
        os.makedirs(self.work_dir, exist_ok=True)

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
    out=$(grep "^$readname" "$sam_path" | awk '{{ if ($9 == 0) print $10 }}')
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

            for step, [sample, chr1, bp1, chr2, bp2] in enumerate(data[head:tail]):
                jun_path = self.jun_dic[sample]
                cmd = cmd_template.format(chr1=chr1, bp1=bp1, chr2=chr2, bp2=bp2,
                                          jun_path=jun_path, out_path=out_path)
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                out, err = p.communicate()
                results.append({'line': head + step + 1, 'ret': p.returncode,
                                'err': err.decode(), 'cnt': out.decode(),
                                'sample': sample, 'chr1': chr1, 'bp1': bp1, 'chr2': chr2, 'bp2': bp2})
            sender.send(results)

        jobs = []
        pipes = []
        for i in range(parallel_num):
            receiver, sender = multiprocessing.Pipe(False)
            mp = multiprocessing.Process(target=task, args=(i, sender))
            jobs.append(mp)
            pipes.append(receiver)
            mp.start()

        # run jobs
        start = time.time()
        for job in jobs:
            job.join()
        self.elapsed_time['preprocess'] = f'{time.time() - start}'

        # results
        data = []
        has_err = False
        for pipe in pipes:
            results = pipe.recv()
            for r in results:
                if r['ret'] != 0 or len(r['err']) != 0 or not r['cnt'].isdecimal():
                    print(f'[Error] {str(r)}')
                    has_err = True
                r.pop('ret')
                r.pop('err')
                data.append(r)
        if has_err:
            exit(1)

        return data

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
        start = time.time()
        out, err = p.communicate()
        if p.returncode != 0:
            print('[Error] at blat runtime')
        self.elapsed_time['blat'] = f'{time.time() - start}'

    def filter(self, data):
        import time
        start = time.time()
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
grep '^>' | {inp_file} |  sed -e 's/^>//'
'''.format(work_dir=self.work_dir, inp_file=self.files['preproc'])
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
        readnames, err = p.communicate()
        if p.returncode != 0:
            print('[Error] at grep runtime')
        #for d in data:
        #    sample = d['sample']
        #    print(sample)
        self.elapsed_time['filter'] = f'{time.time() - start}'

    def run(self):
        data = self.get_data()
        data = self.preprocess(data)
        print(f'Preprocess elapsed time: {float(self.elapsed_time["preprocess"]):.3f} [s]')
        self.cat()
        self.blat()
        print(f'Blat elapsed time: {float(self.elapsed_time["blat"]):.3f} [s]')
        self.filter(data)
        print(f'Filter elapsed time: {float(self.elapsed_time["filter"]):.3f} [s]')
