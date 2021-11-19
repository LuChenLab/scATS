#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Preprocess of UTR bed and BAM, to get the data for model
"""
import gzip
import os
import random
import re
from typing import Dict, List

import pysam

from src.bam import Bam
from src.loci import BED, Reads
from src.logger import log
from src.progress import custom_progress


def load_ats(path: str) -> Dict:
    u"""
    load ats modeling output data
    :param path: the path to ats modeling output file
    :return dict: keys -> utr region; values -> list of splice region
    """

    beds = {}
    header = None
    progress = custom_progress(io=True)
    with progress:
        task_id = progress.add_task(
            f"Reading... ", total=os.path.getsize(path))

        with open(path) as r:
            for _, line in enumerate(r):
                progress.update(task_id, advance=len(str.encode(line)))
                line = line.strip().split("\t")

                if not header:
                    header = line
                    continue

                data = {i: j for i, j in zip(header, line)}

                chrom = data["utr"].split(":")
                chrom = chrom[0]

                strand = "+" if data["utr"].endswith(":+") else "-"
                key = data["gene_name"]

                if "," in key:
                    continue

                alpha = data["infered_sites"].split(",")

                if key not in beds.keys():
                    beds[key] = set()
                try:
                    for x in alpha:
                        if x != "":
                            x = int(float(x))

                            beds[key].add(BED(
                                chrom, x - 1, x, strand,
                                "{}_{}".format(data["gene_name"], data["transcript_name"]), str(
                                    len(beds) + 1)
                            ))
                except Exception as e:
                    log.debug(e)
                    pass
    return beds


def __extract_information__(line: str):
    u"""
    extract information from gff of gtf file
    :param line: string after column 8 of gff|gtf file, eg: ID=xxx;Name=xxx or gene_id "xxx"; gene_name "xxx"
    :return: dict
    """
    data = {}

    for i in line.split(";"):
        i = i.strip()
        if not i:
            continue

        tmp = re.split(r"[= ]", i)

        tmp_info = tmp[1].strip() if ":" not in tmp[1] else tmp[1].split(":")[
            1].strip()
        data[tmp[0]] = tmp_info.replace("\"", "")

    return data


def load_gtf(path: str) -> Dict:
    u"""
    load ats modeling output data
    :param path: the path to ats modeling output file
    :param min_ats: the minimum number of ATS locate in single UTR
    :return dict: keys -> utr region; values -> list of splice region
    """

    beds = {}

    progress = custom_progress(io=True)
    with progress:
        task_id = progress.add_task(
            f"Reading... ", total=os.path.getsize(path))

        with open(path) as r:
            for _, line in enumerate(r):
                progress.update(task_id, advance=len(str.encode(line)))

                if line.startswith('#'):
                    continue

                line = line.strip().split()
                if len(line) <= 8:
                    continue

                if line[2] != "transcript" and not "RNA" in line[2]:
                    continue

                strand = line[6]

                chrom, start_pos, end_pos = line[0], int(line[3]), int(line[4])
                info = __extract_information__(" ".join(line[8:]))
                gene = info.get("gene_name", info.get("Parent"))

                if gene not in beds.keys():
                    beds[gene] = set()

                beds[gene].add(BED(
                    chrom, start_pos, end_pos, strand,
                    name=info.get("transcript_name", info.get("Name")),
                    record_id=""
                ))

    return beds


def load_utr(path: str, debug: bool = False) -> List[BED]:
    u"""
    Load extracted UTR from bed file
    :param path: path to bed file
    :return list of BED objects
    """

    res = []

    progress = custom_progress(io=True)

    with progress:
        task_id = progress.add_task(
            f"Reading...  ", total=os.path.getsize(path))

        with gzip.open(path, "rt") if path.endswith("gz") else open(path) as r:
            for line in r:

                if debug and len(res) > 50:
                    return res
                progress.update(task_id, advance=len(str.encode(line)))

                b = BED.create(line)
                res.append(b)

    if not debug:
        random.seed(42)
        random.shuffle(res)
    return res


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


def load_reads(bam: List[Bam], region: BED, remove_duplicate_umi: bool = False):
    u"""
    Load reads, keys -> R1; values -> R2
    Only both R1

    @2021.06.10 remove dict from this function, use list and sort to increase function speed, with disadvantage like more memory usage

    :params bam: list of bam files
    :params region:
    :return generator: generate r1 and r2
    """
    not_paired = False
    for b in bam:
        umis = set()
        r1s, r2s = [], []

        r = pysam.AlignmentFile(str(b)) if isinstance(b, Bam) else b

        # check the chromosome format
        if region.chromosome not in r.references:
            chromosome = region.chromosome.replace(
                "chr", "") if region.chromosome.startswith("chr") else "chr" + region.chromosome

            if chromosome not in r.references:
                log.warn(
                    f"{region.chromosome} or {chromosome} not present in bam file")
                continue

            region.chromosome = chromosome

        # fetch with until_eof is faster on large bam file according to
        # https://pysam.readthedocs.io/en/latest/faq.html
        for rec in r.fetch(region.chromosome, region.start, region.end, until_eof=True):
            if not rec.is_paired:
                not_paired = True

            if rec.is_unmapped or rec.is_qcfail or rec.mate_is_unmapped:
                continue

            if __get_strand__(rec) != region.strand:
                continue

            # only use the required barcodes for analysis
            if not __is_barcode_exists__(b, rec):
                continue

            # filter duplicate umi
            if remove_duplicate_umi:
                if not rec.has_tag("UB") or rec.get_tag("UB") in umis:
                    continue
                umis.add(rec.get_tag("UB"))

            if rec.is_read1:
                r1s.append(rec)
            else:
                r2s.append(rec)

        if isinstance(b, str):
            r.close()

        if not_paired:
            for r in r1s:
                r = Reads.create(r, single_end=True)
                if r:
                    yield r, None

            for r in r2s:
                r = Reads.create(r, single_end=True)
                if r:
                    yield r, None

            return

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
                r1 = Reads.create(r1)
                r2 = Reads.create(r2)
                if r1 and r2 and region.is_cover(r1, 0) and region.is_cover(r2, 0):
                    yield r1, r2

                i += 1


if __name__ == '__main__':
    pass
