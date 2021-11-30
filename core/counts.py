#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Functions used to count the single cell level expression and calculate PSI
"""
import gzip
from multiprocessing import Process, Queue
from typing import Dict, List

from src.bam import Bam
from src.expression import Expr
from src.progress import custom_progress
from src.reader import load_ats, load_gtf


def count_consumer(bam_files: List[Bam], input_queue: Queue, output_queue: Queue, gtf: bool):
    u"""
    Count ATS
    """
    while True:
        data = input_queue.get()
        utr, regions = data

        res = {}
        for r in regions:
            row_id = f"{utr}_{r.to_str()}" if not gtf else f"{utr}-{r.name}"

            if row_id not in res.keys():
                res[row_id] = {}

            for b in bam_files:
                for cb, ub in b.reads(region = r):
                    col_id = f"{cb}_{b.alias}" if b.alias else cb

                    if col_id not in res[row_id].keys():
                        res[row_id][col_id] = set()
                    
                    res[row_id][col_id].add(ub)
        
        output_queue.put(res)


def count_bulk_consumer(bam_files: List[Bam], input_queue: Queue, output_queue: Queue, gtf: bool):
    u"""
    Count ATS from bulk bam
    """
    while True:
        data = input_queue.get()
        utr, regions = data

        res = {}
        for r in regions:
            row_id = f"{utr}_{r.to_str()}" if not gtf else f"{utr}-{r.name}"

            if row_id not in res.keys():
                res[row_id] = {}

            for b in bam_files:
                res[row_id][b.alias] = b.reads_bulk(r)
        
        output_queue.put(res)


def counts(bams: Dict, output: str, ats: str, processes: int = 1, bulk: bool = False) :
    u"""
    count 
    """
    bam_files = []
    for i, j in bams.items():
        bam_files.append(Bam(i, j))

    gtf = ats.endswith("gtf")

    ats = load_gtf(ats) if gtf else load_ats(ats)

    input_queue = Queue()
    output_queue = Queue()

    # generate consumers
    consumers = []
    for _ in range(processes):
        p = Process(
            target=count_consumer if not bulk else count_bulk_consumer,
            args=(
                bam_files,
                input_queue,
                output_queue,
                gtf,
            )
        )
        p.daemon = True
        p.start()
        consumers.append(p)

    for utr, regions in ats.items():
        input_queue.put([utr, regions])

    progress = custom_progress()

    with progress:
        w = gzip.open(output, "wt+") if output.endswith(".gz") else open(output, "w+")
        
        task = progress.add_task("Counting...", total=len(ats))

        while not progress.finished:
            res = output_queue.get()

            for row, data in res.items():
                for col, val in data.items():
                    if not bulk:
                        val = len(val)
                        
                    if val:
                        w.write(f"{row}\t{col}\t{val}\n")
                    w.flush()
            progress.update(task, advance=1)

        w.close()


def psi(mtx: str, output: str):
    u"""
    Calculate PSI base on count matrix
    :params mtx: count matrix
    """

    # cell -> utr -> sum
    expr = Expr.get_psi(mtx)
        
    w = gzip.open(output, "wt+") if output.endswith(".gz") else open(output, "w+")
    for i in expr:
        w.write("\t".join([str(x) for x in i]) + "\n")
    
    w.close()



if __name__ == '__main__':
    pass
