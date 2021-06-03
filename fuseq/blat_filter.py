import os
import shutil
from fuseq.base import Base
from fuseq.timer import Timer

class BlatFilter(Base):
    def __init__(self, params, breakinfo):
        super().__init__()
        self.params = params
        self.breakinfo = breakinfo

    def __update_chr_bp(self, i, d):
        return d[i]['chr1'], d[i]['chr2'], int(d[i]['bp1']), int(d[i]['bp2'])

    def __write_to_file(self, fws, breakinfo, readname, seq, poses_filt, other_info):
        chr1 = breakinfo['chr1']
        bp1 = breakinfo['bp1']
        strand1 = breakinfo['strand1']
        gene1 = breakinfo['gene1']
        junc1 = breakinfo['junc1']
        chr2 = breakinfo['chr2']
        bp2 = breakinfo['bp2']
        gene2 = breakinfo['gene2']
        junc2 = breakinfo['junc2']
        strand2 = breakinfo['strand2']
        mfline = breakinfo['linenr']
        readname = readname[readname.index('_') + 1:]  # Remove the leading unique number

        w1 = f'{readname} fusionLineNr={mfline}\n'
        w2 = ' '.join((chr1, bp1, strand1, gene1, junc1)) + '\n'
        w3 = ' '.join((chr2, bp2, strand2, gene2, junc2)) + '\n'
        w4 = seq + '\n'
        w = w1 + w2 + w3 + w4
        if poses_filt:
            fw = fws[0]
            itemitems = [[f'{" "*(s-1)}{seq[s-1:e]}', f'[{s},{e}]', f'[{bps},{bpe}]', f'chr{chr}', f'{strand}']
                         for s, e, _, bps, bpe, chr, strand in poses_filt]
            strstrs = [[items[0] for items in itemitems], [items[1] for items in itemitems],
                       [items[2] for items in itemitems], [items[3] for items in itemitems]]
            lenlens = [[len(str) for str in strs] for i, strs in enumerate(strstrs)]
            maxlens = [max(lens) for lens in lenlens]
            spsps = [[maxlens[i] - le for le in lens] for i, lens in enumerate(lenlens)]
            w5 = ''
            for i, items in enumerate(itemitems):
                w5 += f'{items[0]}{" "*spsps[0][i]} {items[1]}{" "*spsps[1][i]} {items[2]}{" "*spsps[2][i]} {items[3]}{" "*spsps[3][i]} {items[4]}\n'
            w += w5
        else:
            fw = fws[1]
        #w += str(other_info) + '\n'
        fw.write(w + '\n')
        #print(w, end='')
        #print('-' * 126)

    def __pos_filter(self, poses):
        '''Filter base on qstat-qend range'''
        if len(poses) < 2:
            return []

        starts = [r[0] for r in poses]
        idx_asc = sorted(range(len(starts)), key=lambda k: starts[k])
        poses_asc = [poses[i] for i in idx_asc]  # poses[i] = (start, end, one_or_two, bp_start, bp_end, chr, strand)

        cur_idx = 0
        max_idx = len(idx_asc) - 1
        results = [0] * len(idx_asc)  # if 0/1, unnecessary/necessary elements
        while cur_idx < max_idx:
            s, e = poses_asc[cur_idx][0], poses_asc[cur_idx][1]
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

        poses_filt = [poses[idx_asc[i]] for i, zero_or_one in enumerate(results) if zero_or_one == 1]
        if len(poses_filt) < 2:
            return []

        return poses_filt

    def pos_sort(self, poses_filt, breakinfo):
        if not poses_filt:
            return poses_filt, breakinfo

        # Check
        one_or_two_list = [p[2] for p in poses_filt]
        if len(set(one_or_two_list)) != 2:
            print('[Error] poses_filt has some problem')
            print(poses_filt)
            exit(1)

        # Sort position and break information
        # if one_or_two == 1, sequence1 is displayed first before sequence2
        # if one_or_two == 2, sequence2 is displayed first before sequence1
        starts = [p[0] for p in poses_filt]
        idx_min = starts.index(min(starts))
        is_one_first = True if poses_filt[idx_min][2] == 1 else False
        if is_one_first:
            idxes = [i for i, one_or_two in enumerate(one_or_two_list) if one_or_two == 1] + \
                    [i for i, one_or_two in enumerate(one_or_two_list) if one_or_two == 2]
        else:
            idxes = [i for i, one_or_two in enumerate(one_or_two_list) if one_or_two == 2] + \
                    [i for i, one_or_two in enumerate(one_or_two_list) if one_or_two == 1]
        new_poses_filt = [poses_filt[idx] for idx in idxes]

        if new_poses_filt[0][2] == 1:
            new_breakinfo = breakinfo
        else:  # Reverse
            b = breakinfo.copy()
            b['chr2'], b['bp2'], b['strand2'], b['chr1'], b['bp1'], b['strand1'] = \
            b['chr1'], b['bp1'], b['strand1'], b['chr2'], b['bp2'], b['strand2']
            new_breakinfo = b

        return new_poses_filt, new_breakinfo

    def __filter_and_write(self, fws, breakinfo, readname, seq, poses, other_info):
        poses_filt = self.__pos_filter(poses)
        poses_filt, breakinfo = self.pos_sort(poses_filt, breakinfo)
        self.__write_to_file(fws, breakinfo, readname, seq, poses_filt, other_info)
        poses.clear()

    # [1] bp1 and bp2 are in the range of Tstart and Tend
    # [2] chr1 and chr2 are the same as Tname
    # [3] No range if only one item satisfies both [1] and [2]
    # [4] In the case of three or more items that satisfy [1] and [2]
    #     when the range of [Qstart2,Qend2] is wider than that of [Qstart1,Qend1], [Qstart2,Qend2] is given priority
    #     [Qstart1,Qend1] is not displayed
    @Timer('filter')
    def run(self):
        # Get read names
        cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
grep '^>' {inp_file} | sed -e 's/^>//'
'''.format(work_dir=self.params.work_dir, inp_file=self.files['coll'])
        readnames = self._run_cmd(cmd, 'grep')
        readnames = readnames.split('\n')

        # Check
        cnts = [int(d.get('cnt')) for d in self.breakinfo]
        assert(len(readnames) == sum(cnts))

        # Create a dictionary with readname in key and breakinfo index in value
        i_rn = 0
        rn2idx = dict()
        for i, cnt in enumerate(cnts):
            for _ in range(cnt):
                rn2idx[readnames[i_rn]] = i
                i_rn += 1

        # Full path
        blat_path = f'{self.params.work_dir}/{self.files["blat"]}'
        coll_path = f'{self.params.work_dir}/{self.files["coll"]}'
        filtmatch_path = f'{self.params.work_dir}/{self.files["filtmatch"]}'
        filtmiss_path = f'{self.params.work_dir}/{self.files["filtmiss"]}'
        fuseq_path = self.params.fuseq_path

        # Open
        fr_collect = open(coll_path, 'r')
        fr_blat = open(blat_path, 'r')
        fw_filtmatch = open(filtmatch_path, 'w')
        fw_filtmiss = open(filtmiss_path, 'w')
        fw_filts = [fw_filtmatch, fw_filtmiss]

        # Filter and write
        is_current = True
        bp_start_extn = self.params.bp_start_extn
        bp_end_extn = self.params.bp_end_extn
        for row in fr_collect:
            # Target read
            cur_readname = row.rstrip('\n')[1:]
            cur_seq = fr_collect.readline().rstrip('\n')
            cur_chr1, cur_chr2, cur_bp1, cur_bp2 = self.__update_chr_bp(rn2idx[cur_readname], self.breakinfo)
            cur_breakinfo = self.breakinfo[rn2idx[cur_readname]]

            poses = []
            other_info = []
            while True:
                if is_current:
                    s = fr_blat.readline().rstrip('\n').split('\t')
                    if s == ['']:  # Check if file pointer is last
                        break
                    strand = s[8]           # strand
                    readname = s[9]         # qname
                    pos_start = int(s[11])  # qstart
                    pos_end = int(s[12])    # qend
                    chr = s[13]             # tname
                    bp_start = int(s[15])   # tstart
                    bp_end = int(s[16])     # tend
                    pos_start_adj = pos_start + 1
                    bp_start_adj = bp_start - bp_start_extn
                    bp_end_adj = bp_end + bp_end_extn  # NOTE: +1 increases the number of matches to Genomon results
                if readname == cur_readname:
                    # Filter based on chr and tstart-tend range
                    if chr == cur_chr1 and bp_start_adj <= cur_bp1 <= bp_end_adj:
                        poses.append((pos_start_adj, pos_end, 1, bp_start + 1, bp_end, chr, strand))
                    elif chr == cur_chr2 and bp_start_adj <= cur_bp2 <= bp_end_adj:
                        poses.append((pos_start_adj, pos_end, 2, bp_start + 1, bp_end, chr, strand))
                    if chr == cur_chr1:
                        other_info.append((bp_start_adj, bp_end_adj, 1))
                    elif chr == cur_chr2:
                        other_info.append((bp_start_adj, bp_end_adj, 2))
                    is_current = True
                else:
                    # Write one read
                    self.__filter_and_write(fw_filts, cur_breakinfo, cur_readname, cur_seq, poses, other_info)
                    is_current = False
                    break
        if is_current:
            # Write the remaining read
            self.__filter_and_write(fw_filts, cur_breakinfo, cur_readname, cur_seq, poses, other_info)

        # Close
        fr_collect.close()
        fr_blat.close()
        fw_filtmatch.close()
        fw_filtmiss.close()

        # Concat filtmatch and filtmiss
        filtmatch_exists = os.path.exists(filtmatch_path)
        filtmiss_exists = os.path.exists(filtmiss_path)
        if filtmatch_exists and filtmiss_exists:
            cmd = '''\
#!/bin/bash
set -eu
cd {work_dir}
cat {filtmatch} {filtmiss} > {fuseq}
'''.format(work_dir=self.params.work_dir, fuseq=fuseq_path, filtmatch=filtmatch_path, filtmiss=filtmiss_path)
            self._run_cmd(cmd, 'cat_filtmatch_filtmiss')
        elif filtmatch_exists:
            shutil.move(filtmatch_path, fuseq_path)
        elif filtmiss_exists:
            shutil.move(filtmiss_path, fuseq_path)
