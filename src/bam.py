#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.11.16
"""
import gzip
import os
import sys

import pysam

from src.loci import BED
from src.logger import log


class Bam(object):

    def __init__(self, path: str, alias=None):
        self.path = os.path.abspath(os.path.realpath(path))
        self.alias = os.path.basename(path) if not alias else alias

        if not self.check_bam(self.path):
            log.error(f"{self.path} is no a valid bam")
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

        log.debug(f"Load barcodes from {path}")
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

    def reads(self, region: BED, cell_tag: str = "CB", umi_tag: str = "UB"):
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

    def reads_bulk(self, region: BED):
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
    pass
