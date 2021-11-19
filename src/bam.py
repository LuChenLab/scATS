#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.11.16
"""
import gzip
import os
import sys

import pysam

# from src.loci import BED
# from src.logger import log


class Bam(object):

    def __init__(self, path: str, alias=None):
        self.path = os.path.abspath(os.path.realpath(path))
        self.alias = os.path.basename(path) if not alias else alias

        if not self.check_bam(self.path):
            # log.error(f"{self.path} is no a valid bam")
            sys.exit(1)

        self.barcodes = self.load_barcodes(os.path.join(os.path.dirname(
            self.path), "filtered_feature_bc_matrix/barcodes.tsv.gz"))

    def __str__(self):
        return self.path

    def __hash__(self):
        return hash(self.path)

    @staticmethod
    def check_bam(path: str) -> bool:
        try:
            with pysam.AlignmentFile(path, require_index=True) as r:
                pass
        except IOError as err:
            log.info(f"try to create index for {path}")

            try:
                pysam.index(path)
            except Exception as err:
                log.error(err)
                sys.exit(err)
        except Exception as err:
            return False
        return True

    @staticmethod
    def load_barcodes(path: str) -> dict:
        u"""
        load required barcodes

        :params path
        """
        res = {}

        if not os.path.exists(path):
            return res

        r = gzip.open(path, "rt") if path.endswith(".gz") else open(path)

        for line in r:
            line = line.strip()

            key = line[:min(3, len(line))]
            temp = res.get(key, set())
            temp.add(line)
            res[key] = temp

        r.close()

        return res

    @property
    def has_barcode(self) -> bool:
        return len(self.barcodes) > 0

    def check(self, barcode: str):
        if not self.has_barcode:
            return True

        key = barcode[:min(3, len(barcode))]
        barcodes = self.barcodes.get(key, set())
        return barcode in barcodes

    def reads(self, region, cell_tag: str = "CB", umi_tag: str = "UB"):
        u"""
        generator to get all cell barcode and umi barcode from specific region
        :params region: target region
        :params cell_tag: the cell barcode tag in 10x bam
        :params umi_tag: the umi tag in 10x bam
        :return cell barcode and umi
        """
        with pysam.AlignmentFile(self.path) as r:
            for record in r.fetch(region.chromosome, region.start, region.end):
                if record.is_qcfail or record.is_unmapped:
                    continue

                if record.has_tag("NH") and record.get_tag("NH") > 1:
                    continue

                if record.has_tag(cell_tag) and record.has_tag(umi_tag):
                    yield record.get_tag(cell_tag), record.get_tag(umi_tag)

    def reads_bulk(self, region):
        u"""
        generator to get all cell barcode and umi barcode from specific region
        :params region: target region
        :return the nubmer of reads
        """
        count = 0
        with pysam.AlignmentFile(self.path) as r:
            for record in r.fetch(region.chromosome, region.start, region.end):
                if record.is_qcfail or record.is_unmapped:
                    continue

                if record.has_tag("NH") and record.get_tag("NH") > 1:
                    continue

                count += 1
        return count


if __name__ == "__main__":
    import pysam

    bam = Bam("/mnt/raid64/ATS/alignments/cellranger/hs/vdj_v1_hs_nsclc_5gex/outs/possorted_genome_bam.bam")


    def __get_strand__(read: pysam.AlignedSegment) -> str:
        u"""
        determine the reads strand

        :params read: 
        """

        if read.is_paired:
            if read.is_read1 and read.is_reverse:
                return "-"
            elif read.is_read2 and not read.is_reverse:
                return "-"
            return "+"

        return "-" if read.is_reverse else "+"


    def __is_barcode_exists__(b: Bam, rec: pysam.AlignedSegment) -> bool:
        u"""
        check whether this read contains required barcode
        :params barcodes: a collection of required barcodes
        :params rec: 
        """
        if not b.has_barcode:
            return True

        if not rec.has_tag("CB"):
            return False

        cb = rec.get_tag("CB")
        return b.check(cb)

    r1s, r2s = [], []
    with pysam.AlignmentFile(str(bam)) as r:
        for rec in r.fetch("1", 1074057, 1074557, until_eof=True):
           
            if rec.is_unmapped or rec.is_qcfail or rec.mate_is_unmapped:
                continue

            if __get_strand__(rec) != "-":
                continue

            # only use the required barcodes for analysis
            if not __is_barcode_exists__(bam, rec):
                continue

            if rec.is_read1:
                r1s.append(rec)
            else:
                r2s.append(rec)

    r1s = sorted(r1s, key=lambda x: x.query_name)
    r2s = sorted(r2s, key=lambda x: x.query_name)

    i, j = 0, 0
    while i < len(r1s) and j < len(r2s):
        r1 = r1s[i]
        r2 = r2s[j]

        if r1.query_name < r2.query_name:
            i += 1
        elif r1.query_name > r2.query_name:
            j += 1
        else:

            print(r1.query_name, r2.query_name)

            i += 1


    pass
