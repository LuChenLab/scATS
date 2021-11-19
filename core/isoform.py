#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.05.06 by Zhang

This script contains the function to infer isoforms
"""
import gzip
import os
import re
from functools import lru_cache, reduce
from typing import List

import numpy as np
import pysam
from src.convert import Coordinate, TreeNode, Window, WinList
from src.loci import BED, GTF, GenomicLoci
from src.logger import log
from src.progress import custom_progress

np.seterr(divide='ignore', invalid='ignore')

#################### CODE STRUCTURE #######################
# Part 1: Isoform mapping, tree construction
# Part 2: isoform processing and assignment
# Part 3: ATS mixture

################### Part 1: Isoform mapping #####################################


def map_iso_to_gene(iso_list: WinList, gene_list: WinList):
    """
    :param iso_list: exon windows of the isoform on gene
    :param gene_list: smallest windows across isoforms by intersecting all exons of all isoforms on the gene
    :return: a WinList (same size as gene_list), with indexes relative to the isoform (local index)
    """
    if not iso_list.is_sorted:
        iso_list.sort()

    if not gene_list.is_sorted:
        gene_list.sort()
    assert iso_list.is_sorted, "iso list is not sorted"
    assert gene_list.is_sorted, "gene list is not sorted"
    res_win_list = WinList()
    curr_exon_j = 0  # jth exon on isoform
    local_iso_index = 0
    for win in gene_list:
        if curr_exon_j >= len(iso_list):
            res_win_list.append(Window())
            continue
        if win <= iso_list[curr_exon_j]:
            res_win_list.append(win.shift_start(local_iso_index))
            local_iso_index += len(win)
            if iso_list[curr_exon_j].end == win.end:
                curr_exon_j += 1
        else:
            res_win_list.append(Window())
    return res_win_list


def build_tree(isowins_list: List[WinList]):
    """
    :param isowins_list: each element is an aligned WinList of the same window across different isoforms
    :return: a tree
    """

    def dfs(i_left, i_right):
        root = TreeNode()
        if i_left == i_right - 1:
            root.winlist = alignwins_list[i_left]
            root.left = None
            root.right = None
        else:
            mid = (i_left + i_right) // 2
            tmplist = [WinList(isowins[i_left:i_right]).get_range()
                       for isowins in isowins_list]
            root.winlist = WinList(tmplist)
            root.left = dfs(i_left, mid)
            root.right = dfs(mid, i_right)
        return root

    # each element is an aligned WinList of the same window across different isoforms
    alignwins_list = [WinList(wins) for wins in zip(*isowins_list)]

    return dfs(0, len(alignwins_list))


def concatenate_winlist(winlist1: WinList, winlist2: WinList):
    res_list = [concatenate_windows(w1, w2)
                for w1, w2 in zip(winlist1, winlist2)]
    return WinList(res_list)


def concatenate_windows(win_left: Window, win_right: Window) -> Window:
    if win_left and win_right and win_left.end == win_right.start:
        return Window(win_left.start, win_right.end)
    else:
        return Window()


def query(root: TreeNode, query_winlist: WinList, iso_index: int):
    """
    map the query_winlist to isoform, get the relative index of the query window on all isoforms
    :param root: root of the tree
    :param query_winlist: winlist of a read mapped to the "iso_index"th isoform
    :param iso_index: index of given isoform, used for locating the window
    :return: mapped windows on all isoforms
    """

    @lru_cache(maxsize=1024)
    def query_window(root: TreeNode, qwin: Window, i_iso: int):
        if root is None or root.winlist[i_iso].is_empty():
            return none_winlist
        if not qwin <= root.winlist[i_iso]:
            return none_winlist
        if root.left is None and root.right is None:  # leaf case
            local_st = qwin.start - root.winlist[i_iso].start
            local_win = Window(local_st, local_st + len(qwin))
            res_winlist = WinList()
            for win in root.winlist:
                if not win:
                    res_winlist.append(Window())
                else:
                    res_winlist.append(
                        local_win.shift_start(win.start + local_st))
            return res_winlist
        # root is inner node, one child might be intron
        if root.left.winlist[i_iso].end != Window.DEL_VAL:
            mid = root.left.winlist[i_iso].end
        else:
            mid = root.right.winlist[i_iso].start

        if qwin.end <= mid:
            return query_window(root.left, qwin, i_iso)
        elif qwin.start >= mid:
            return query_window(root.right, qwin, i_iso)
        else:
            qwin_left = Window(qwin.start, mid)
            qwin_right = Window(mid, qwin.end)
            left_winlist = query_window(root.left, qwin_left, i_iso)
            right_winlist = query_window(root.right, qwin_right, i_iso)

            return WinList([concatenate_windows(lw, rw) for lw, rw in zip(left_winlist, right_winlist)])

    n_isoform = len(root.winlist)
    none_winlist = WinList([Window() for _ in range(n_isoform)])

    res_list = [query_window(root, win, iso_index) for win in query_winlist]
    res_winlist = reduce(concatenate_winlist, res_list)

    return res_winlist


# map a window to gene
def map_to_gene(root: TreeNode, query_winlist: WinList, iso_index: int):
    """
    map the query_winlist to gene, get the relative index of the query window on the gene
    :param root: root of the tree, the first exon must be the full gene
    :param query_winlist: a WinList on the "iso_index"th isoform
    :param iso_index: index of given isoform, used for locating the window
    :return: mapped windows on the gene
    """

    @lru_cache(maxsize=1024)
    def map_window(root: TreeNode, qwin: Window, i_iso: int):
        if root is None or root.winlist[i_iso].is_empty():
            return none_list
        if not qwin <= root.winlist[i_iso]:
            return none_list
        if root.left is None and root.right is None:  # leaf case
            local_st = qwin.start - root.winlist[i_iso].start
            local_win = Window(local_st, local_st + len(qwin))
            gene_win = root.winlist[0]
            res_winlist = WinList(
                [local_win.shift_start(gene_win.start + local_st), ])
            return res_winlist
        # root is inner node, one child might be intron
        if root.left.winlist[i_iso].end != Window.DEL_VAL:
            mid = root.left.winlist[i_iso].end
        else:
            mid = root.right.winlist[i_iso].start

        if qwin.end <= mid:
            return map_window(root.left, qwin, i_iso)
        elif qwin.start >= mid:
            return map_window(root.right, qwin, i_iso)
        else:
            qwin_left = Window(qwin.start, mid)
            qwin_right = Window(mid, qwin.end)
            left_list = map_window(root.left, qwin_left, i_iso)
            right_list = map_window(root.right, qwin_right, i_iso)
            return left_list+right_list

    none_list = WinList()

    res_list = [map_window(root, win, iso_index) for win in query_winlist]
    res_winlist = reduce(lambda a, b: a + b, res_list)

    return res_winlist


################### Part 2: isoform processing and assignment #####################################
# modify exons of isoform such that ATS falls into an exon
def proc_isoform(ats_pos: int, iso_wins: WinList):
    """
    :param iso_wins_list: list of isoforms, 0th element is gene range
    :return: a list of isoforms with exons modified
    """
    # step 1: check if [max(0, ats-50), ats+50] falls in exon
    # step 1.1: if yes => no special treatment
    # step 1.2: if no  => merge [max(0, ats-50), min(ats+200, exon_start)] with the nearest first exon
    ats_win = Window(ats_pos, ats_pos + 1)
    i = 0
    while i < len(iso_wins) and iso_wins[i] << ats_win:
        i += 1
    if i > 0:  # ATS is after the first exon
        return None
    first_exon_ind = i  # first downstream exon

    ats_win = Window(max(0, ats_pos - 50), ats_pos + 50)
    if ats_win <= iso_wins[first_exon_ind]:
        return iso_wins

    expanded_ats_win = Window(
        max(0, ats_pos - 50), min(ats_pos + 200, iso_wins[first_exon_ind].start + 1))

    # expanded ats window must overlap with isoform
    if expanded_ats_win << iso_wins[first_exon_ind]:
        return None

    iso_wins.append(expanded_ats_win)
    iso_wins.is_sorted = False
    return iso_wins.merge()  # sorted wins


# calculate relative positions of reads on each isoform
def get_read_relative_pos(iso_wins_list: List[WinList], r1_wins_list: List[WinList], r2_wins_list: List[WinList],
                          read_labels):
    # break all isoforms into atomic non-overlapping windows
    all_wins = WinList([win for iso_wins in iso_wins_list for win in iso_wins])
    all_win_on_gene = all_wins.split()

    # map all isoforms to the gene (atomic non-overlapping windows)
    atom_isowins_list = []
    for iso_wins in iso_wins_list:
        atom_isowins_list.append(map_iso_to_gene(iso_wins, all_win_on_gene))

    # Modified by Zhang at 2021.05.08
    # some cases the iso_wins cannot map to gene
    if len(atom_isowins_list) < 1:
        return None, None

    # build the tree
    tree = build_tree(atom_isowins_list)

    # get relative position of R1 and R2 on the gene
    n_frag = len(read_labels)
    n_iso = len(iso_wins_list)
    read_st_mat = np.zeros((n_frag, n_iso), dtype=int)
    read_en_mat = np.zeros((n_frag, n_iso), dtype=int)

    for i in range(n_frag):
        r1_wins = r1_wins_list[i]
        r2_wins = r2_wins_list[i]
        res1 = query(tree, r1_wins, read_labels[i])
        res2 = query(tree, r2_wins, read_labels[i])
        for j in range(n_iso):
            if res1[j] and res2[j]:
                read_st_mat[i, j] = res1[j].start
                read_en_mat[i, j] = res2[j].end

    return read_st_mat, read_en_mat


# given an ATS position, calculate the posterior distribution of different isoforms
# if return all zeros, it hints this may be a novel ATS
def assign_isoform(ats_pos: int, iso_wins_list: List[WinList], r1_wins_list: List[WinList], r2_wins_list: List[WinList],
                   read_labels, mu_frag=400, sd_frag=30, min_frag_len=150):
    """
    :param ats_pos: position of given ATS
    :param iso_wins_list: a list of isoform, each is a WinList of Windows (exons).
            the first must be the full gene Window(0, utr_length)
    :param r1_wins_list: r1 winlist, a consecutive read can map to multiple windows on the gene
    :param r2_wins_list: r2 winlist
    :param read_labels: isoform index of each pair-read, 0 is the full gene
    :param mu_frag: mean of fragment length, e.g. 400
    :param sd_frag: sd of fragment length, e.g. 30
    :param min_frag_len: minimum fragment length, fragments shorter than this are equally assigned to all isoforms
    :return: posterior probability of different isoforms, e.g. np.array([0, 0.3, 0.7, 0])
            the first exon will always have weight 0
    """
    # extend ATS and merge with exons

    # 0th element of iso_wins_list is the whole gene window i.e. [Window(0,utr_length),]
    # return None if ats_pos is not near the first exon
    iso_wins_list = [proc_isoform(ats_pos, iso_wins)
                     for iso_wins in iso_wins_list]
    n_orig_iso = len(iso_wins_list)
    valid_iso_inds = np.array([i for i, w in enumerate(iso_wins_list) if w])

    # in case ats_pos is not close to the start position of any isoform, i.e. ats_pos only compatible with the whole gene
    if len(valid_iso_inds) == 1:
        return np.zeros(n_orig_iso, dtype='float')

    # remove empty isoforms
    iso_wins_list = [w for i, w in enumerate(iso_wins_list) if w]

    # change read labels from all isoforms to selected isoforms
    label_map = {k: i for i, k in enumerate(valid_iso_inds)}
    # unselected isoforms will be mapped to the first exon
    read_labels = [label_map.get(k, 0) for k in read_labels]

    # get relative positions of the reads on each isoform
    read_st_mat, read_en_mat = get_read_relative_pos(
        iso_wins_list, r1_wins_list, r2_wins_list, read_labels)

    # Modified by Zhang at 2021.05.08
    # check whether failed to get relative positions
    if read_st_mat is None and read_en_mat is None:
        return np.zeros(n_orig_iso, dtype='float')

    # assign each pair-end read to different isoforms
    n_frag, n_iso = read_st_mat.shape
    read_len_mat = read_en_mat - read_st_mat
    post_prob_mat = np.zeros((n_frag, n_iso), dtype='float')
    # all-zero rows in read_len_mat are properly handled
    post_prob_mat[:, 1:] = cal_post_prob(mu_frag, sd_frag, read_len_mat[:, 1:])

    # special treatment to short fragments, overwrite previous calculation
    # if the fragment is less then the minimum fragment size (150bp),
    # they should be assigned to all valid isoforms equally
    # min_frag_len = 150
    # only need to be treated if mapped to multiple isoforms
    flag_arr1 = np.sum(read_len_mat[:, 1:], axis=1) > 1
    # fragment length less than min_frag_len
    flag_arr2 = np.max(read_len_mat[:, 1:], axis=1) < min_frag_len
    short_frag_inds = np.logical_and(flag_arr1, flag_arr2)
    tmpmat = read_len_mat[short_frag_inds, 1:] > 0
    post_prob_mat[short_frag_inds, 1:] = tmpmat / \
        np.sum(tmpmat, axis=1)[:, np.newaxis]

    iso_post_prob = np.sum(post_prob_mat, axis=0)

    iso_post_prob = iso_post_prob / np.sum(iso_post_prob)  # 对应于该ATS的各isoform比例

    res_post_prob = np.zeros(n_orig_iso, dtype='float')
    res_post_prob[valid_iso_inds] = iso_post_prob

    return res_post_prob


# calculate the posterior probability of each DNA fragment on different isoforms
def cal_post_prob(mu_frag_size: int, sd_frag_size: int, len_iso_mat: np.ndarray):
    """
    :param mu_frag_size: mean of fragment size
    :param sd_frag_size: std of fragment size
    :param len_iso_mat: n_frag x n_iso numpy matrix, each element is the length of DNA fragment on each isoform
    :return: n_frag x n_iso numpy matrix, each element is the posterior of a given fragment on an isoform
    """

    def norm_pdf(x, mean, sd):
        prob_density = (1 / (np.sqrt(2 * np.pi) * sd)) * \
            np.exp(-0.5 * ((x - mean) / sd) ** 2)
        return prob_density

    max_arr = np.max(len_iso_mat, axis=1)
    val_inds = max_arr > 0

    prob_mat = np.zeros(len_iso_mat.shape, dtype="float")
    prob_mat[val_inds, :] = norm_pdf(
        len_iso_mat[val_inds, :], mu_frag_size, sd_frag_size)
    # len=0 means reads cannot be mapped to isoform
    prob_mat[len_iso_mat == 0] = 0

    prob_mat[val_inds, :] = prob_mat[val_inds, :] / \
        np.sum(prob_mat[val_inds, :], 1)[:, np.newaxis]

    return prob_mat


class GTFUtils(object):
    u"""
    object to handle the GTF related process

    including 
    - gtf format validation
    - gtf file indexing
    - load gtf record
    - convert absolute gtf to relative
    """

    def __init__(self, path: str, iterator: bool = True):
        u"""
        init the obj with path to gtf file
        """
        assert os.path.exists(path), f"{path} not exists"
        assert self.is_gtf(path), f"{path} is not gtf"
        self.gtf = os.path.abspath(path)

        if iterator:
            self.gtf = self.index_gtf()

    @staticmethod
    def is_gtf(infile):
        u"""
        check if input file is gtf
        :param infile: path to input file
        :return: status code, 0 -> invalid gtf, 1 -> valid gtf file 
        """
        if infile is None:
            return False

        is_gtf = 0
        try:
            r = gzip.open(infile, "rt") if infile.endswith(
                ".gz") else open(infile)
        except TypeError as err:
            log.error("failed to open %s", infile)
            exit(err)

        for line in r:
            if line.startswith("#"):
                continue

            lines = re.split(r"\s+", line)

            if len(lines) < 8:
                break

            if re.search(
                r"([\w-]+ \"[\w.\s\-%,:]+\";? ?)+",
                " ".join(lines[8:])
            ):
                is_gtf += 1

            break

        r.close()

        return is_gtf

    def index_gtf(self, sort_gtf=True, retry=0):
        u"""
        Created by ygidtu
        Extract only exon tags and keep it clean
        :param sort_gtf: Boolean value, whether to sort gtf file first
        :param retry: only try to sort gtf once
        :return path to compressed and indexed bgzipped gtf file
        """
        if self.gtf is None:
            return None

        gtf = self.is_gtf(self.gtf)

        if gtf % 10 != 1:
            raise ValueError(
                "gtf file required, %s seems not a valid gtf file" % self.gtf)

        index = False
        if gtf // 10 > 0:
            output_gtf = self.gtf
        else:
            output_gtf = self.gtf + ".gz"
        if not os.path.exists(output_gtf) or not os.path.exists(output_gtf + ".tbi"):
            index = True

        elif os.path.getctime(output_gtf) < os.path.getctime(output_gtf) or \
                os.path.getctime(output_gtf) < os.path.getctime(output_gtf):
            index = True

        # 2018.12.21 used to handle gtf not sorted error
        if sort_gtf and retry > 1:
            raise OSError(
                "Create index for %s failed, and trying to sort it failed too" % self.gtf)
        elif sort_gtf:
            data = []

            # log.info("Sorting %s" % self.gtf)

            old_gtf = self.gtf
            self.gtf = re.sub("\.gtf$", "", self.gtf) + ".sorted.gtf"

            output_gtf = self.gtf + ".gz"

            if os.path.exists(self.gtf) and os.path.exists(output_gtf):
                return output_gtf

            try:
                w = open(self.gtf, "w+")
            except IOError as err:
                w = open("/tmp/sorted.gtf")

            with open(old_gtf) as r:
                for line in r:
                    if line.startswith("#"):
                        w.write(line)
                        continue

                    lines = line.split()

                    if len(lines) < 1:
                        continue

                    data.append(
                        GenomicLoci(
                            chromosome=lines[0],
                            start=lines[3],
                            end=lines[4],
                            strand=lines[6],
                            gtf_line=line
                        )
                    )

            for i in sorted(data):
                w.write(i.gtf_line)

            w.close()

        if index:
            # log.info("Create index for %s", self.gtf)
            try:
                pysam.tabix_index(
                    self.gtf,
                    preset="gff",
                    force=True,
                    keep_original=True
                )
            except OSError as err:

                if re.search("could not open", str(err)):
                    raise err

                # log.error(err)
                # log.error("Guess gtf needs to be sorted")
                return self.index_gtf(sort_gtf=True, retry=retry + 1)

        return output_gtf

    def read_transcripts(self, region: BED) -> Coordinate:
        u"""
        Read transcripts from tabix indexed gtf files
        The original function check if the junctions corresponding to any exists exons, I disable this here
        :param self.gtf: path to bgzip gtf files (with tabix index), only ordered exons in this gtf file
        :param region: splice region
        :param retry: if the gtf chromosome and input chromosome does not match. eg: chr9:1-100:+ <-> 9:1-100:+
        :param genome: path to genome fasta file
        :return: SpliceRegion and gene region
        """

        res = {}
        try:
            # log.info("Reading from %s" % self.gtf)

            with pysam.Tabixfile(self.gtf, 'r') as gtf_tabix:
                relevant_exons_iterator = gtf_tabix.fetch(
                    region.chromosome,
                    region.start - 1,
                    region.end + 1,
                    parser=pysam.asGTF()
                )

                # min_exon_start, max_exon_end, exons_list = float("inf"), float("-inf"),  []
                transcripts = []
                gene_ids = set()
                gene_records = {}
                for line in relevant_exons_iterator:
                    try:
                        if line.feature == "transcript" and line.strand == region.strand:
                            transcripts.append(BED(
                                line.contig,
                                line.start,
                                line.end,
                                line.strand,
                                line.transcript_id,
                                line.gene_id
                            ))
                            gene_ids.add(line.gene_id)

                        if line.feature == "gene" and line.strand == region.strand:
                            gene_records[line.gene_id] = BED(
                                line.contig,
                                line.start,
                                line.end,
                                line.strand,
                                "",
                                line.gene_id
                            )

                    except IndexError as err:
                        log.error(err)
                        # print(err)

                # if there are multiple genes at same utr, just give up
                if len(gene_ids) > 1:
                    return None

                # 2nd round iteration by transcripts to retrive all exons of transcript
                for t in transcripts:
                    exons = []
                    for line in gtf_tabix.fetch(t.chromosome, t.start - 1, t.end + 1, parser=pysam.asGTF()):
                        try:
                            if line.feature == "exon" and line.transcript_id == t.name:
                                exons.append(BED(
                                    line.contig,
                                    line.start,
                                    line.end,
                                    line.strand,
                                    line.exon_id,
                                    line.transcript_id
                                ))

                        except IndexError as err:
                            log.error(err)
                            # print(err)

                    res[t] = exons
        except ValueError as err:
            log.warn(err)
            return None
        return Coordinate(gene=gene_records[list(gene_ids)[0]], isoforms=res)
    
    def read_gtf(self, debug: bool = False):
        path = self.gtf
        # read gtf
        progress = custom_progress(io = True)
        plus, minus = [], []
        gene = None
        transcripts = {}

        with progress:
            task_id = progress.add_task(f"Reading...", total=os.path.getsize(path))

            with gzip.open(path, "rt") if path.endswith("gz") else open(path, "r") as r:
                for line in r:
                    progress.update(task_id, advance=len(str.encode(line)))

                    if line.startswith("#"):
                        continue

                    rec = GTF.create(line)

                    if rec.source == "gene":
                        if gene != rec:
                            if gene is not None:
                                if gene.strand == "+":
                                    plus.append(Coordinate(gene, transcripts))
                                else:
                                    minus.append(Coordinate(gene, transcripts))

                            gene = rec
                            transcripts = {}

                    if rec.source == "exon":
                        if rec.transcript_name not in transcripts.keys():
                            transcripts[rec.transcript_name] = []

                        transcripts[rec.transcript_name].append(rec)
                    
                    if debug and (len(plus) + len(minus)) > 30:
                        break

        if gene is not None:
            if gene.strand == "+":
                plus.append(Coordinate(gene, transcripts))
            else:
                minus.append(Coordinate(gene, transcripts))
            
        del gene, transcripts

        # merge gene regions for further filtering
        res = []
        for temp in [plus, minus]:
            temp = sorted(temp, key=lambda x: [x.chromosome, x.start, x.end])
            
            if not temp:
                continue

            curr = temp[0]
            for i in temp[1:]:
                if curr & i:
                    curr += i
                else:
                    res.append(curr)
                    curr = i
            res.append(curr)

        return sorted(res, key=lambda x: [x.chromosome, x.start, x.end])


if __name__ == '__main__':
    pass
