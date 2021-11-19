#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.05.07 by Zhang

This script contains the objects and functions related to the coordinate convertion
"""

from typing import Dict, List, Optional, Tuple

import pysam

from src.loci import BED, Reads, Region
from src.reader import load_reads


class Window:
    DEL_VAL = 1 << 31

    def __init__(self, start: int = DEL_VAL, end: int = DEL_VAL, is_read1: bool = False):
        """
        :type start: int
        :type end: int
        """
        self.start = start
        self.end = end
        self.is_read1 = is_read1

    def __key(self) -> Tuple[int, int]:
        return self.start, self.end

    def __hash__(self):
        return hash(self.__key())

    def is_empty(self):
        return self.start == self.end

    def __bool__(self):
        return self.start != Window.DEL_VAL

    def __len__(self):
        return self.end - self.start

    def __str__(self):
        return f'{self.start}  {self.end}  {len(self)}'

    #################### window vs window -> bool #####################
    def __lshift__(self, other):
        return self.end <= other.start

    def __lt__(self, other):
        return self.start < other.start < self.end < other.end

    def __le__(self, other):
        return other.start <= self.start and self.end <= other.end

    def __eq__(self, other):
        if isinstance(other, Window):
            return self.start == other.start and self.end == other.end
        else:
            return False

    def __ge__(self, other):
        return self.start <= other.start and other.end <= self.end

    def __gt__(self, other):
        return other.start < self.start < other.end < self.end

    def __rshift__(self, other):
        return other.end <= self.start

    # overlap
    def __ne__(self, other):
        return not (other.end <= self.start or self.end <= other.start)

    # adjacent
    def adj(self, other):
        return self.end == other.start or self.start == other.end

    #################### window vs point operation #####################
    # if a point falls within the window
    def __contains__(self, point):
        if self.is_empty():
            return False
        return self.start <= point < self.end

    def shift_start(self, start):
        if self.start == Window.DEL_VAL:
            return self
        else:
            return Window(start, start + self.end - self.start)

    @classmethod
    def create(cls, region: Reads):
        u"""
        convert retion to Window
        """
        return cls(region.start, region.end, is_read1=region.is_read1)


class WinList(list):
    def __init__(self, *args, is_sorted=False):
        super().__init__(*args)
        self.is_sorted = is_sorted

    def __str__(self):
        return "\n".join(['start  end  length'] + [str(w) for w in self.__iter__()])

    def append(self, __object: Window) -> None:
        super().append(__object)
        self.is_sorted = False

    def rm_empty_win(self):
        return WinList([w for w in self if w], is_sorted=self.is_sorted)

    def sort(self) -> None:
        if not self.is_sorted:
            super().sort(key=lambda win: (win.start, win.end))
            self.is_sorted = True

    def rm_duplicate(self):
        res = WinList(set(self))
        res.sort()
        return res

    # merge overlapping windows
    def merge(self):
        winlist = self.rm_empty_win()
        assert len(winlist) > 0
        if not winlist.is_sorted:
            winlist.sort()
        res_list = WinList([winlist[0]])
        for win in winlist:
            # note that windows are sorted first by start, then by end
            curr_win = res_list[-1]
            if curr_win.start <= win.start <= curr_win.end:
                res_list[-1] = Window(curr_win.start, win.end)
            else:
                res_list.append(win)
        res_list.is_sorted = True
        return res_list

    # split the winList into tiny non-overlapping intervals, including regions that are not covered
    # e.g. [(0,3), (2,4), (6,8)] => [(0,2), (2,3), (3,4), (4,6), (6,8)]
    def split(self):
        winlist = self.rm_empty_win()
        if len(winlist) == 0:
            return winlist
        if not winlist.is_sorted:
            winlist.sort()
        boarder = set()
        for win in winlist:
            if not win:
                continue
            boarder.add(win.start)
            boarder.add(win.end)
        boarder_arr = sorted(list(boarder))
        winlist = [Window(i, j) for i, j in zip(boarder_arr, boarder_arr[1:])]
        return WinList(winlist, is_sorted=True)

    # return the left most position of the windows list
    # mainly for building trees, note that there may be Empty windows in this list, non-empty windows are sorted
    def get_left(self):
        for win in self:
            if win.start != Window.DEL_VAL:
                return win.start
        return Window.DEL_VAL

    # return the right most position of the windows list
    def get_right(self):
        for win in reversed(self):
            if win.end != Window.DEL_VAL:
                return win.end
        return Window.DEL_VAL

    def get_range(self):
        return Window(self.get_left(), self.get_right())


class TreeNode:
    def __init__(self, left=None, right=None, winlist: WinList = WinList()):
        """
        :type left: TreeNode
        :type right: TreeNode
        :type winlist: a list of windows, each
        """
        self.left = left
        self.right = right
        self.winlist = winlist  # ith window corresponds to range of this tree node


class Coordinate(object):
    u"""
    object to handle the coordanite convertion of isoforms from reference file
    """

    __slots__ = ('gene', 'isoforms', 'ids', 'bams')

    def __init__(self, gene: Region, isoforms:  Dict):
        u"""
        init
        :params gene: gene region
        :params isoforms: transcript id -> list of exons
        """

        self.gene = gene
        self.isoforms = isoforms
        self.ids = self.__generate_isoform_idx__(isoforms)
        self.bams = []

    @property
    def chromosome(self):
        return self.gene.chromosome
    
    @property
    def start(self):
        return self.gene.start

    @property
    def end(self):
        return self.gene.end
    
    @property
    def strand(self):
        return self.gene.strand

    def __and__(self, other):
        return self.gene & other.gene

    def __add__(self, other):
        gene = self.gene + other.gene

        isoforms = self.isoforms
        isoforms.update(other.isoforms)

        return Coordinate(gene, isoforms)

    def set_bams(self, bams: List[str]):
        u"""
        check input bam files
        """
        for bam in bams:
            try:
                with pysam.AlignmentFile(str(bam)) as r:
                    pass
            except Exception as err:
                raise FileNotFoundError(
                    f"{bam} is not a valid bam file: {err}")

        self.bams = bams

    @staticmethod
    def __generate_isoform_idx__(isoforms: Dict):
        u"""
        assign each isoform a unique index id
        """
        data = []
        for iso, exons in isoforms.items():
            data.append([iso, exons])

        data = sorted(data, key=lambda x: [x[1][0].start, x[1][-1].end])

        return [x[0] for x in data]

    @property
    def relative(self) -> List:
        u"""
        convert absolute coordinate to relative
        :return [[id, start, end], [id, start, end]]
        """

        res = [[0, 0, len(self.gene)]]

        for idx, iso in enumerate(self.ids):
            idx += 1
            exons = self.isoforms[iso]

            for e in exons:
                res.append([idx, e.start - self.gene.start,
                           e.end - self.gene.start])

        return res

    @property
    def winList(self) -> List[WinList]:
        relative = self.relative
        return [WinList([Window(x[1], x[2]) for x in relative if x[0] == i]) for i in range(len(self.ids) + 1)]

    @property
    def is_duplicate(self):
        return self.gene.is_duplicate
    
    def get(self, idx: int) -> Optional[str]:
        u"""
        get transcript id by index
        """
        if idx == 0:
            return self.gene.id

        if idx > len(self.ids):
            return None

        return self.ids[idx - 1].name

    def reads(self, region: BED, remove_duplicate_umi: bool = False) -> Dict:
        u"""
        load reads and convert and assign isoforms
        """
        return load_reads(
            self.bams, region,
            remove_duplicate_umi=remove_duplicate_umi
        )

    def __is_loci_inside_of_gene__(self, pos: int) -> bool:
        u"""
        Check whether a position is located inside of gene region
        """
        return self.gene.start <= pos < self.gene.end

    def __relative_in_gene__(self, pos: int) -> int:
        u"""
        Convert a position to relative to gene
        """
        return pos - self.gene.start

    def reads_to_relative(self, reads: Reads, return_winlist: bool = False):
        u"""
        convert reads coord to relative
        """
        res = []
        if not reads:
            return res
        exons = reads.exon_parts

        for i in range(0, len(exons), 2):
            start = exons[i]
            end = exons[i + 1]

            if self.__is_loci_inside_of_gene__(start) and self.__is_loci_inside_of_gene__(end):
                start = self.__relative_in_gene__(start)
                end = self.__relative_in_gene__(end)
                res.append(Reads(
                    ref=self.gene.chromosome,
                    start=start, end=end,
                    strand=self.gene.strand,
                    is_read1=reads.is_read1,
                    name=reads.name
                ))
            elif self.__is_loci_inside_of_gene__(start):
                res.append(Reads(
                    ref=self.gene.chromosome,
                    start=self.__relative_in_gene__(start),
                    end=self.__relative_in_gene__(start) + end - start,
                    strand=self.gene.strand,
                    is_read1=reads.is_read1,
                    name=reads.name
                ))
            elif self.__is_loci_inside_of_gene__(end):
                res.append(Reads(
                    ref=self.gene.chromosome,
                    start=self.__relative_in_gene__(end) - end + start,
                    end=self.__relative_in_gene__(end),
                    strand=self.gene.strand,
                    is_read1=reads.is_read1,
                    name=reads.name
                ))

        if return_winlist:
            return WinList([Window.create(x) for x in res])
        return res

    def assign(self, reads: Reads, utr_length: int = 500):
        u"""
        wether all reads parts locate in exons of same transcript
        :param reads: list of Reads
        :param paired
        """
        assigned = []

        # convert reads to relative coord and cast into Region
        reads = sorted(self.reads_to_relative(reads))

        if not reads:
            return assigned

        # cast transcripts to Region
        transcripts = []
        for idx, i in enumerate(self.relative):
            if i[0] > 0:
                st, en = i[1], i[2]

                if idx == 0:
                    st = max(st - utr_length // 2, 1)
                
                if idx == len(self.relative) - 1:
                    en = en + utr_length // 2

                transcripts.append(BED(
                    self.gene.chromosome,
                    st, en,
                    self.gene.strand,
                    record_id=i[0],
                    name=i[0]
                ))

        transcripts = sorted(transcripts)

        # make sure the reads and transcripts are sorted
        matches = {}
        i, j = 0, 0
        while i < len(reads) and j < len(transcripts):
            current_reads = reads[i]
            current_trans = transcripts[j]

            if current_reads < current_trans:
                i += 1
            elif current_reads > current_trans:
                j += 1
            else:
                if current_trans.is_cover(current_reads, 5):
                    matches[current_trans.id] = 1 + \
                        matches.get(current_trans.id, 0)
                i += 1

        for idx, match in matches.items():
            if match == len(reads):
                assigned.append(idx)

        return assigned

    def utr_per_transcript(self, span: int = 500):
        utrs = []
        for iso, exons in self.isoforms.items():
            exons = sorted(exons, key=lambda x: [x.chromosome, x.start, x.end])

            site = exons[0].start if self.gene.strand == "+" else exons[-1].end

            utrs.append(BED(
                self.gene.chromosome,
                max(site - span, 1),
                site + span,
                self.gene.strand,
                name=iso,
                record_id=self.gene.gene_name
            ))

        res = []

        if utrs:
            utrs = sorted(utrs, key=lambda x: [x.chromosome, x.start, x.end])
            curr_utr = utrs[0]

            for i in utrs[1:]:
                if curr_utr & i:
                    curr_utr =  curr_utr + i
                else:
                    res.append(curr_utr)
                    curr_utr = i
            res.append(curr_utr)
        return res


if __name__ == '__main__':
    from rich import print
    gene = BED("1", 1000, 4000, "+", "", "")

    isos = {
        "iso1": [
            BED("1", 1000, 1150, "+", "", ""),
            BED("1", 1116, 2100, "+", "", ""),
            BED("1", 1823, 3100, "+", "", "")
        ],
        "iso2": [
            BED("1", 1102, 3020, "+", "", ""),
            BED("1", 1203, 3300, "+", "", ""),
            BED("1", 2604, 3560, "+", "", "")
        ],
        "iso3": [
            BED("1", 3405, 4000, "+", "", ""),
            BED("1", 3706, 4000, "+", "", ""),
            BED("1", 3843, 4000, "+", "", "")
        ]
    }

    rel = ReferenceIsoform(gene, isos)

    print(rel.ids)
    print(rel.relative)
    print(rel.get(4))
