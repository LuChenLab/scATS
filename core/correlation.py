#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by Zhang at 2021.06.16

Core function to calculate co-expression
"""
import os
import re

from typing import Optional

from scipy.stats import pearsonr, spearmanr

from src.expression import Expr
from src.progress import custom_progress


class Barcodes(object):
    u"""
    Used this to colllect the barcodes with group info and check the existance
    """

    def __init__(self, path: str):
        self.data = {}

        progress = custom_progress(io = True)

        idx = -1
        with progress:
            task_id = progress.add_task("Load barcodes...", total = os.path.getsize(path))
            with open(path) as r:
                for line in r:
                    progress.update(task_id, advance=len(line.encode()))
                    line = line.strip().split()

                    for i, j in enumerate(line):
                        if re.search(r"[ATCG]+(-1)?", j):
                            idx = i
                            break

                    if idx < 0:
                        raise ValueError("is there any barcode in %s ?" % line)

                    key = line[idx + 1]
                    barcode = line[idx]

                    if key not in self.data.keys():
                        self.data[key] = set()
                    
                    self.data[key].add(barcode)

    def exists(self, group: str, barcode: str) -> bool:
        if not group in self.data.keys():
            return False

        prefix = barcode[:min(3, len(barcode))]
        return barcode in self.data[group].get(prefix, set())

    def __iter__(self):
        for key, val in self.data.items():
            yield key, val


def __correction_between_ats__(
    first_ats: dict, 
    second_ats: dict, 
    distance: int,
    expr: Expr,
    corr_func, 
    barcodes = None
):
    u"""
    as namy says

    :params first_ats: list of row_id
    :params second_ats: list of row_id
    :params distance: the maxmimum distance between two UTR
    :params expr: all the expression data
    :params barcodes: all the barcodes
    :params corr_func: 
    :params a set of barcodes
    """
    atslist1 = sorted(first_ats.keys())
    atslist2 = sorted(second_ats.keys())

    i, j = 0, 0 # backup index, current index
    res = []
    for ats1 in atslist1:

        for j in range(i, len(atslist2)):
            ats2 = atslist2[j]

            if ats1.start - ats2.end > distance:
                i = -1
                break
            elif ats2.start - ats1.end > distance:
                continue
            else:
                if i < 0:
                    i = j
                ats1_expr = expr.get_expr(first_ats[ats1], barcodes = barcodes)
                ats2_expr = expr.get_expr(second_ats[ats2], barcodes = barcodes)

                if len(ats1_expr) < 2 and len(ats2_expr) < 2:
                    continue

                r, p = corr_func(ats1_expr, ats2_expr)

                res.append(f"{ats1}\t{ats2}\t{r}\t{p}")

    return res


def corr(mtx: str, output: str, group_info: Optional[str] = None, distance: int = 1000, pearson: bool = True, barcode: str = None):
    u"""
    calculate correction of ATS count

    :params mtx: path to count file
    :params output: path to output file
    :params distance: the maxmimum distance between two UTR
    :params pearson: whether to calculate pearson or spearman correction
    :params barcode: path to list of barcodes to use
    """

    corr_func = pearsonr if pearson else spearmanr

    if group_info:
        group_info = Barcodes(group_info)

    expr = Expr.create(mtx, barcode)
    expr.sort()

    progress = custom_progress()

    with progress:
        task_id = progress.add_task("Computing... ", total=len(expr))

        with open(output, "w+") as w:
            
            if group_info:
                for group, barcodes in group_info:
                    for i in range(len(expr)):
                        last_utr = expr.utrs[i]

                        for j in range(i+1, len(expr.utrs)):
                            curr_utr = expr.utrs[j]

                            # if utr is too far away, skip
                            if curr_utr.chromosome != last_utr.chromosome or curr_utr.start - last_utr.end > distance:
                                break

                            if last_utr.strand != curr_utr.strand:
                                for r in __correction_between_ats__(
                                    first_ats=expr.get_row_ids(last_utr),
                                    second_ats=expr.get_row_ids(curr_utr),
                                    distance=distance,
                                    expr=expr,
                                    corr_func=corr_func,
                                    barcodes = barcodes
                                ):
                                    w.write(group + "\t" + r + "\n")
            else:

                for i in range(len(expr)):
                    last_utr = expr.utrs[i]

                    for j in range(i+1, len(expr.utrs)):
                        curr_utr = expr.utrs[j]

                        # if utr is too far away, skip
                        if curr_utr.chromosome != last_utr.chromosome or curr_utr.start - last_utr.end > distance:
                            break

                        if last_utr.strand != curr_utr.strand:
                            for r in __correction_between_ats__(
                                first_ats=expr.get_row_ids(last_utr),
                                second_ats=expr.get_row_ids(curr_utr),
                                distance=distance,
                                expr=expr,
                                corr_func=corr_func
                            ):
                                w.write(r + "\n")

                progress.update(task_id, advance=1)


if __name__ == '__main__':
    pass