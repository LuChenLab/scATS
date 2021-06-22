#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Functions used to extract utr from gtf file
"""
import gzip
import os

from rich.progress import track
from src.loci import GTF
from src.progress import custom_progress


def process(gtf: str, output: str, span: int = 500):
    u"""
    extract utr from gtf file
    :param gtf: path to gtf file
    :param output: path to output file
    :param span: the radius of UTR
    :return None
    """
    output = os.path.abspath(output)
    os.makedirs(os.path.dirname(output), exist_ok = True)

    transcripts = {}

    # read gtf
    progress = custom_progress(io = True)
    with progress:
        task_id = progress.add_task(f"Processing...", total=os.path.getsize(gtf))

        with gzip.open(gtf, "rt") if gtf.endswith("gz") else open(gtf, "r") as r:
            for line in r:
                progress.update(task_id, advance=len(str.encode(line)))

                if line.startswith("#"):
                    continue

                rec = GTF.create(line)

                if rec.source == "exon":
                    if rec.transcript_id not in transcripts.keys():
                        transcripts[rec.transcript_id] = []

                    transcripts[rec.transcript_id].append(rec)

    # get first exons
    first_exons = []
    for exons in track(transcripts.values(), description=f"Collecting first exons"):
        exons = sorted(exons)
        exon = exons[0]
        site = exon.start
        if exons[0].strand == "-":
            exon = exons[-1]
            site = exon.end
        
        exon = GTF(
            chromosome=exon.chromosome,
            start=max(site - span, 1),
            end=site + span,
            strand=exon.strand,
            attrs=exon.attrs,
            source = exon.source
        )
        first_exons.append(exon)
    first_exons = sorted(first_exons, key=lambda x: (x.chromosome, x.start, x.end))

    # merging
    with gzip.open(output, "wt+") if output.endswith("gz") else open(output, "w+") as w:
        curr_exon = first_exons[0]
        for i in track(range(1, len(first_exons)), description=f"Writing {os.path.basename(output)}"):
            if curr_exon & first_exons[i]:
                curr_exon = curr_exon + first_exons[i]
            else:
                w.write(f"{curr_exon.bed}\n") 
                curr_exon = first_exons[i]
        w.write(f"{curr_exon.bed}\n")
                        

if __name__ == '__main__':
    pass
