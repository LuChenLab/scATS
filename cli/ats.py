#!/usr/bin/envpython3
# -*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Contians all the parameters and command line params handler
"""
from multiprocessing import Process, Queue, cpu_count
from typing import List, Optional

import click

from src.logger import init_logger, log
from src.progress import custom_progress
from src.reader import check_bam, load_barcodes
from core.isoform import GTFUtils

from core.model import AtsModel, Parameters


class ATSParams(object):
    u"""
    This a a param handler for ATS inferrence
    """

    def __init__(
        self,
        gtf: str,
        bam: List[str],
        n_max_ats: int = 5,
        n_min_ats: int = 1,
        min_ws: float = 0.01,
        max_unif_ws: float = 0.1,
        max_beta: int = 50,
        fixed_inference_flag: bool = False,
        debug: bool = False,
        remove_duplicate_umi: bool = False,
        min_reads: int = 10,
        utr_length: int = 500
    ):
        u"""
        init this function

        :params bam: path to bam files or list of bams
        """
        log.info("Load UTR")

        self.index, self.bam = None, None

        gtf = GTFUtils(gtf, iterator = False)
        self.gtf = gtf.read_gtf()
        self.bam = self.check_path(bam)

        self.barcodes = {b: load_barcodes(b.replace(".bam", ".barcode")) for b in self.bam}

        if not self.bam:
            raise ValueError("Please input valid bam files")

        self.n_max_ats = n_max_ats
        self.n_min_ats = n_min_ats
        self.min_ws = min_ws
        self.max_unif_ws = max_unif_ws
        self.max_beta = max_beta
        self.fixed_inference_flag = fixed_inference_flag
        self.debug = debug
        self.remove_duplicate_umi = remove_duplicate_umi
        self.min_reads = min_reads
        self.utr_length = utr_length

    @staticmethod
    def check_path(bams: List[str]) -> List[str]:
        u"""
        check input bam files
        """
        for bam in bams:
            if not check_bam(bam):
                log.error(f"{bam} is not a valid bam file")
                exit(1)
        return bams

    def __str__(self) -> str:
        res = []
        for i in dir(self):
            if i.startswith("__"):
                continue

            if "function" in str(type(getattr(self, i))):
                continue

            res.append(f"- {i}: {getattr(self, i)}")

        return "\n".join(res)

    def __len__(self):
        return len(self.gtf)

    @staticmethod
    def __format_reads_to_relative__(reads, utr) -> List[int]:
        u"""
        format list of reads into start site array
        """
        st_arr = []

        utr_site = utr.start if utr.strand == "+" else utr.end

        # Only iter reads1
        for r, _ in reads:
            if utr.start <= r.start <= r.end <= utr.end:
                site = r.start if utr.strand == "+" else r.end
                st_arr.append(site - utr_site if utr.strand ==
                              "+" else utr_site - site)

        return st_arr

    def keys(self) -> List[str]:
        res = ["utr", "reference_id", "reference_name", "number_of_reads", "infered_sites"]
        res += Parameters.keys()
        return res

    def format_res(self, utr, res: Parameters, number_of_reads: int = 0) -> str:
        u"""
        as name says format ATS model results to meaningful str
        """
        site = utr.start if utr.strand == "+" else utr.end

        sites = [str(site + x if utr.strand == "+" else site - x)
                 for x in res.alpha_arr]

        data = [
            f"{utr.chromosome}:{utr.start}-{utr.end}:{utr.strand}",
            utr.id,
            utr.name,
            str(number_of_reads),
            ",".join(sites),
            res.to_res_str()
        ]
        return "\t".join(data)

    def run(self, idx: int) -> Optional[str]:
        u"""
        Factory function to execute the ATS model and format results
        """
        gene = self.gtf[idx]
        gene.set_bams(self.bam)

        utrs = []
        reads = {}
        for utr in gene.utr_per_transcript(span = self.utr_length // 2):
            utrs.append(utr)

            for r1, r2 in gene.reads(utr, remove_duplicate_umi = self.remove_duplicate_umi):
                assign = set(gene.assign(r1)) & set(gene.assign(r2))
                
                for a in assign:
                    if a not in reads:
                        reads[a] = []

                    reads[a].append([r1, r2])

        for idx, rds in reads.items():    
            if len(rds) < self.min_reads:
                continue
            utr = utrs[idx - 1]
            st_arr = self.__format_reads_to_relative__(rds, utr)

            if len(st_arr) <= 1:
                continue

            m = AtsModel(
                n_max_ats=self.n_max_ats,
                n_min_ats=self.n_min_ats,
                st_arr=st_arr,
                utr_length=len(utr),
                min_ws=self.min_ws,
                max_unif_ws=self.max_unif_ws,
                max_beta=self.max_beta,
                fixed_inference_flag=self.fixed_inference_flag
            )

            try:
                res = m.run()
                if res:
                    yield f"{self.format_res(utr, res, len(m.st_arr))}" 
            except Exception as err:
                if self.debug:
                    log.exception(err)
                    exit(0)
                else:
                    log.error(err)


def consumer(input_queue: Queue, output_queue: Queue, params: ATSParams):
    u"""
    Multiprocessing consumer to perform the ATS core function
    :param input_queue: multiprocessing.Queue to get the index
    :param output_queue: multiprocessing.Queue to return the results
    :param params: the parameters for ATS model
    """

    while True:
        gene = input_queue.get()

        if gene is None:
            break

        for res in params.run(gene):
            output_queue.put(res)


@click.command()
@click.option(
    "-g", "--gtf",
    type=click.Path(exists=True),
    required=True,
    help=""" The path to gtf file. """
)
@click.option(
    "-o", "--output",
    type=click.Path(),
    required=True,
    help=""" The path to output file. """
)
@click.option(
    "-u", "--utr-length",
    type=int,
    default=500,
    help=""" The length of UTR. """
)
@click.option(
    "--n-max-ats",
    type=click.IntRange(1, 999),
    default=5,
    help=""" The maximum number of ATSs in same UTR. """
)
@click.option(
    "--n-min-ats",
    type=click.IntRange(1, 998),
    default=1,
    help=""" The minimum number of ATSs in same UTR. """
)
@click.option(
    "--min-ws",
    type=float,
    default=0.01,
    help=""" The minimum weight of ATSs. """
)
@click.option(
    "--min-reads",
    type=int,
    default=10,
    help=""" The minimum number of reads in UTR. """
)
@click.option(
    "--max-unif-ws",
    type=float,
    default=0.1,
    help=""" The maximum weight of uniform component. """
)
@click.option(
    "--max-beta",
    type=int,
    default=50,
    help=""" The maximum std for ATSs. """
)
@click.option(
    "--fixed-inference",
    is_flag=True,
    default=True,
    type=click.BOOL,
    help=""" Inference with fixed parameters. """
)
@click.option(
    "-d", "--debug",
    is_flag=True,
    type=click.BOOL,
    help=""" Enable debug mode to get more debugging information, Never used this while running. """
)
@click.option(
    "-p", "--processes",
    type=click.IntRange(1, cpu_count()),
    default=1,
    help=""" How many cpu to use. """
)
@click.option(
    "--remove-duplicate-umi",
    is_flag=True,
    type=click.BOOL,
    help=""" Only kept reads with different UMIs for ATS inference. """
)
@click.argument("bams", nargs=-1, type=click.Path(exists=True), required=True)
def ats(
    gtf: str,
    output: str,
    n_max_ats: int,
    n_min_ats: int,
    utr_length: int,
    min_ws: float,
    min_reads: int,
    max_unif_ws: float,
    max_beta: int,
    fixed_inference: bool,
    processes: int,
    debug: bool,
    remove_duplicate_umi: bool,
    bams: List[str],
):
    u"""
    Inference
    \f

    :param debug: enable debug mode
    """

    init_logger("DEBUG" if debug else "INFO")
    log.info("ATS inference")

    params = ATSParams(
        gtf=gtf,
        bam=bams,
        n_max_ats=n_max_ats,
        n_min_ats=n_min_ats,
        min_ws=min_ws,
        max_unif_ws=max_unif_ws,
        max_beta=max_beta,
        fixed_inference_flag=fixed_inference,
        debug=debug,
        remove_duplicate_umi = remove_duplicate_umi,
        utr_length=utr_length,
        min_reads=min_reads,
    )

    # init queues
    input_queue = Queue()
    output_queue = Queue()

    # generate consumers
    consumers = []
    for _ in range(processes):
        p = Process(
            target=consumer,
            args=(input_queue, output_queue, params,)
        )
        p.daemon = True
        p.start()
        consumers.append(p)


    # producer to assign task
    for i in range(len(params)):
        input_queue.put(i)

    progress = custom_progress()
    with progress:
        task = progress.add_task("Computing...", total=len(params))
        with open(output, "w+") as w:
            header = '\t'.join(params.keys())
            w.write(f"{header}\n")

            while not progress.finished:
                res = output_queue.get(block=True, timeout=None)
                if res:
                    w.write(f"{res}\n")
                    w.flush()

                progress.update(task, advance=1)

    log.info("DONE")


if __name__ == '__main__':
    ats()
