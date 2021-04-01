#import pandas as pd

class Preprocess():
    def __init__(self, genomon):
        self.path_mf = genomon.path_mf  # merge_fusionfusion
        self.path_jun = genomon.path_jun  # junction
        self._proc()

    def _proc(self):
        chrs = list(map(lambda x: str(x), range(23))) + list('XY')

        for _, pa_mf in self.path_mf.items():
            with open(pa_mf, 'r') as f_mf:
                # Skip header line
                f_mf.readline()

                for row_mf in f_mf:

                    # Extract chrs and bps from merge_fusionfusion file
                    row_mf = row_mf.rstrip('\n').split('\t')
                    if row_mf[1] not in chrs or row_mf[4] not in chrs:
                        continue
                    [sample, chr1, bp1, strand1, chr2, bp2, strand2] = row_mf[0:7]
                    bp1 = str(int(bp1) + 1) if strand1 == '+' else str(int(bp1) - 1)
                    bp2 = str(int(bp2) + 1) if strand2 == '+' else str(int(bp2) - 1)

                    # Extract readnames from junction file
                    path_jun = self.path_jun[sample]
                    readnames = []
                    with open(path_jun, 'r') as f_jun:
                        match = False
                        for row_jun in f_jun:
                            r = row_jun.rstrip('\n').split('\t')
                            if (r[0] == chr1 and r[1] == bp1 and r[3] == chr2 and r[4] == bp2) or \
                               (r[3] == chr1 and r[4] == bp1 and r[0] == chr2 and r[1] == bp2):
                                readnames.append(r[9])
                                match = True
                        if not match:
                            print(f'[Error] No data for {chr1} {bp1} {chr2} {bp2} in {path_jun}')
                            exit(1)

                    # Create data for blat input
                    path_sam = path_jun[:-8] + 'sam'
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

            #df = pd.read_csv(pa_mf, sep='\t', header=0,
            #                 usecols=range(7),
            #                 names=['id', 'chr1', 'bp1', 'strand1',
            #                              'chr2', 'bp2', 'strand2'])
            #df = df.astype({'id': str, 'chr1': str, 'bp1': int, 'strand1': str,
            #                           'chr2': str, 'bp2': int, 'strand2': str})
            #df = df.query(f'(chr1 in {str(chrs)}) | (chr2 in {str(chrs)})')

            #d = df.to_()
            #print(d)

            #for idx, row in df.iterrows():
