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
            fw = fws['filtmatch']
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
            fw = fws['filtmiss']
        fw.write(w + '\n')

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

    def __check_pos(self, poses_filt, fws, breakinfo, readname, seq, poses, other_info):
        if not poses_filt:
            return poses_filt

        # Check for error or warning
        is_err = False
        is_war = False
        one_or_two_set = {p[2] for p in poses_filt}
        if len(one_or_two_set) != 2:
            is_err = True
        else:
            one_or_two_set = set()
            for i, pos in enumerate(poses_filt):
                pos_diff = pos[1] - pos[0]
                bp_diff = pos[4] - pos[3]
                if pos_diff == bp_diff:
                    one_or_two_set.add(pos[2])
            if len(one_or_two_set) != 2:
                is_war = True

        # Return if everything is okay
        if not is_err and not is_war:
            return poses_filt

        # Here is error or warning

        # Information about positions
        poses = [p for p in poses_filt if p[2] == 1] + [p for p in poses_filt if p[2] == 2]
        w_sent = []
        for p in poses:
            w1 = f'Type:{p[2]} '
            w2 = f'Qstart:{p[0]} Qend:{p[1]} Tstart:{p[4]} Tend:{p[3]} '
            w3 = f'Qdiff:{p[1]-p[0]} Tdiff:{p[4]-p[3]}'
            w_sent.append(w1 + w2 + w3)
        w_split = [w.split(' ') for w in w_sent]
        w_len = [[len(w) for w in ws] for ws in w_split]
        r_axis0 = range(len(w_split))
        r_axis1 = range(len(w_split[0]))
        w_max = [max([w_len[j][i] for j in r_axis0]) for i in r_axis1]
        w_space = [[w_max[i] - len(w) for i, w in enumerate(ws)] for ws in w_split]
        w_sent = [' '.join([w_split[i][j] + ' ' * w_space[i][j] for j in r_axis1]) for i in r_axis0]
        w = '\n'.join([w.rstrip(' ') for w in w_sent])

        # Print
        if self.params.print_filt_err and (is_err or is_war):
            if is_err:
                print(f'[Warning] Only type{list(one_or_two_set)[0]} exists')
            else:
                print('[Warning] Qpos diff and Tpos diff do not match')
            print(w)

        # Othre information
        w += '\n'
        w += f'readname   | {str(readname)}\n'
        w += f'sequence   | {str(seq)}\n'
        w += f'poses      | {str(poses)}\n'
        w += f'poses_filt | {str(poses_filt)}\n'
        w += f'other_info | {str(other_info)}\n'
        w += '\n'

        # Write to file
        key = 'filterr' if is_err else'filtwar'
        fws[key].write(w)

        return [] if is_err else poses_filt

    def __pos_sort(self, poses_filt, breakinfo):
        if not poses_filt:
            return poses_filt, breakinfo

        # Sort position and break information
        # if one_or_two == 1, sequence1 is displayed first before sequence2
        # if one_or_two == 2, sequence2 is displayed first before sequence1
        one_or_two_list = [p[2] for p in poses_filt]
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
        poses_filt = self.__check_pos(poses_filt, fws, breakinfo, readname, seq, poses, other_info)
        poses_filt, breakinfo = self.__pos_sort(poses_filt, breakinfo)
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
        filterr_path = f'{self.params.work_dir}/{self.files["filterr"]}'
        filtwar_path = f'{self.params.work_dir}/{self.files["filtwar"]}'
        fuseq_path = self.params.fuseq_path

        # Open
        fr_collect = open(coll_path, 'r')
        fr_blat = open(blat_path, 'r')
        fw_filtmatch = open(filtmatch_path, 'w')
        fw_filtmiss = open(filtmiss_path, 'w')
        fw_filterr = open(filterr_path, 'w')
        fw_filtwar = open(filtwar_path, 'w')
        fw_filts = {'filtmatch': fw_filtmatch, 'filtmiss': fw_filtmiss,
                    'filterr': fw_filterr, 'filtwar': fw_filtwar}

        # Filter and write
        is_current = True
        start_extn = self.params.bp_start_extn
        end_extn = self.params.bp_end_extn
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
                    pos_start_plus1 = pos_start + 1  # NOTE: pos_start+1(base1) matches blat result on web
                    bp_start_plus1 = bp_start + 1    # NOTE: bp_start+1(base1) matches blat result on web
                    bp_start_extn = bp_start_plus1 - start_extn
                    bp_end_extn = bp_end + end_extn  # NOTE: end_extn=1 matches Genomon result
                if readname == cur_readname:
                    pos_intvl_ok = pos_end - pos_start == bp_end - bp_start if self.params.check_pos_intvl else True
                    if pos_intvl_ok:
                        # Filter based on chr and tstart-tend range
                        if chr == cur_chr1 and bp_start_extn <= cur_bp1 <= bp_end_extn:
                            poses.append((pos_start_plus1, pos_end, 1, bp_start_plus1, bp_end, chr, strand))
                        elif chr == cur_chr2 and bp_start_extn <= cur_bp2 <= bp_end_extn:
                            poses.append((pos_start_plus1, pos_end, 2, bp_start_plus1, bp_end, chr, strand))
                    if chr == cur_chr1:
                        other_info.append((bp_start_extn, bp_end_extn, 1, cur_bp1))
                    elif chr == cur_chr2:
                        other_info.append((bp_start_extn, bp_end_extn, 2, cur_bp2))
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
        fw_filterr.close()
        fw_filtwar.close()

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
