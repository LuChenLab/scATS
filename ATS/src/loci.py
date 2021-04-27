#!/usr/bin/env python3
#-*- coding:utf-8 -*-
u"""
Created at 2021.04.25 by Zhang

Dedifned genomic loci object
"""

from typing import Optional


class Region(object):
    u"""
    Created at 2021.04.27 by Zhang

    Class handle the genomic regions
    """

    def __init__(self, chromosome: str, start: int, end: int, strand: str):
        assert strand in ("+", "-"), "strand must be + or -"
        assert end >= start, f"end must bigger than start, current -> start: {start}, end: {end}"

        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}:{self.strand}"

    def __hash__(self):
        return hash(self.__str__())

    def __gt__(self, other) -> bool:
        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome
        
        if self.start != self.start:
            return self.start > other.start

        if self.end != self.end:
            return self.end > other.end
        
        return self.strand > other.strand

    def __lt__(self, other) -> bool:
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome
        
        if self.start != self.start:
            return self.start < other.start

        if self.end != self.end:
            return self.end < other.end
        
        return self.strand < other.strand

    def __and__(self, other) -> bool:
        u"""
        override & for overlap checking
        """
        if self.chromosome == other.chromosome and self.strand == other.strand:
            return self.start < other.end and self.end > other.start
        return False

    def __add__(self, other) -> Optional:
        u"""
        override + for merging regions

        :return Region or None: None means there is no any overlap, Region is the merged region
        """
        if self & other:
            return Region(
                chromosome = self.chromosome,
                strand = self.strand,
                start = min(self.start, other.start),
                end = max(self.end, other.end)
            )


class GTF(Region):
    u"""
    Created at 2021.04.27 by Zhang

    Class handle the records in BED format
    """

    __id_label__ = ["gene_id", "gene_name", "transcript_id", "transcript_name", "exon_id"]

    def __init__(
        self, 
        chromosome: str, 
        start: int, end: int, 
        strand: int, 
        source: str,
        attrs: dict = None
    ):
        super(GTF, self).__init__(chromosome, start, end, strand)
        self.attrs = attrs
        self.source = source

        self.ids = {}
        for i in self.__id_label__:
            self.ids[i] = [self.attrs[i]] if i in self.attrs.keys() else []

    def __hash__(self):
        return hash((self.chromosome, self.start, self.end, self.strand, self.source))

    def __add__(self, other):
        if self & other:
            gtf = GTF(
                chromosome=self.chromosome,
                start=min(self.start, other.start),
                end=max(self.end, other.end),
                strand=self.strand,
                attrs = self.attrs,
                source = self.source
            )

            for i in self.__id_label__:
                self.ids[i] += other.ids[i]
                self.ids[i] = list(set(self.ids[i]))
            return gtf

    @classmethod
    def decode_attrs(cls, attrs: str) -> dict:
        u"""
        Decode attrs from gtf files

        :param attrs: gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1";
        :return: {gene_id: "ENSG00000223972"; gene_version: "5", gene_name: "DDX11L1"}
        """
        res = {}

        for i in attrs.split(";"):
            i = i.strip().replace('"', "").replace(";", "")
            values = i.split(" ")

            if len(values) > 1:
                res[values[0].strip()] = values[1].strip()

        return res

    @classmethod
    def create(cls, record: str):
        u"""
        Create GTF object from gtf record
        :param record: 1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; 
        """

        record = record.strip().split("\t")

        source = record[2]
        if "RNA" in record[2]:
            source = "transcript"

        return cls(
            chromosome=record[0], 
            start=int(record[3]), 
            end=int(record[4]),
            strand=record[6],
            attrs = cls.decode_attrs(record[8]) if len(record) >= 8 else {},
            source=source
        )

    @property
    def transcript_id(self) -> str:
        return ",".join(sorted(self.ids["transcript_id"]))

    @property
    def gene_id(self) -> str:
        return ",".join(sorted(self.ids["gene_id"]))

    @property
    def exon_id(self) -> str:
        return ",".join(sorted(self.ids["exon_id"]))

    @property
    def transcript_name(self) -> str:
        return ",".join(sorted(self.ids["transcript_name"]))

    @property
    def gene_name(self) -> str:
        return ",".join(sorted(self.ids["gene_name"]))

    @property
    def bed(self) -> str:
        r_id = ",".join([x for x in [self.gene_id, self.transcript_id, self.exon_id] if x])
        r_name = ",".join([x for x in [self.gene_name, self.transcript_name] if x])

        return f"{self.chromosome}\t{self.start}\t{self.end}\t{r_id}\t{r_name}\t{self.strand}"


class BED(Region):
    u"""
    Created at 2021.04.27 by Zhang

    Class handle the records in BED format
    """

    def __init__(
        self, 
        chromosome: str, 
        start: str, 
        end: str, 
        strand: str, 
        name: str, 
        record_id: str
    ):
        super(BED, self).__init__(chromosome, start, end, strand)
        self.name = name
        self.id = record_id

    def __str__(self) -> str:
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.id}\t{self.name}\t{self.strand}"

    @classmethod
    def create(cls, bed: str):
        u"""
        Create Region obj from bed record

        :param bed: the str of bed record
        """
        bed = bed.strip().split("\t")

        if len(bed) not in [3, 4, 6]:
            raise TypeError(f"Invalid columns, 3, 4 or 6 is required, current: {bed}")

        strand = "+"
        name = "."
        record_id = "."
        
        if len(bed) == 4:
            strand = bed[3]

        if len(bed) == 6:
            record_id = bed[3]
            name = bed[4]
            strand = bed[5]

        return cls(
            chromosome=bed[0],
            start=int(bed[1]),
            end=int(bed[2]),
            record_id=record_id,
            name=name,
            strand=strand
        )
        

class Reads(Region):
    u"""
    Created at 2021.04.27 by Zhang

    Class handle the reads and it's junctions
    """

    def __init__(
        self, 
        ref: str, 
        start: int, 
        end: int, 
        name: str,
        strand: str,
        is_read1: bool, 
        cigar = None
    ):
        u"""
        init the reads object

        :param ref: reference id, chr1, chr2 etc.
        :param start: left-most start
        :param end: right-most site
        :param name: reads name
        :param strand: strand
        :param is_read1: as name says
        :param cigar: cigar tuples, detailes see pysam.AlignmentSegment
        :param paired: SE -> None, PE -> pointer to paired reads 
        """

        super(Region, self).__init__(ref, start, end, strand)

        self.name = name
        self.is_read1 = 

        self.paired = None

    def set_paired(self, v):
        u"""
        setter of paired
        :param v: another Reads object
        """

        self.paired = v
        v.paired = self

    @classmethod
    def create(cls, record: pysam.AlignmentSegment, skip_qc: bool = False):
        u"""
        Create Reads obj from pysam.AlignmentSegment

        :param record: as type
        :param skip_qc: skip QC filtering
        :return if qc is enabled and the records failed qc, then return None
        """
        if record.is_unmapped or record.is_qcfail or not reord.is_proper_pair:
            return None

        if record.has_tag("NH") and record.get_tag("NH") > 1:
            return None

        return cls(
            ref=record.reference_id,
            start=record.reference_start + 1,
            end = record.reference_end + 1,
            strand = cls.__determine_strand__(record),
            name = query_name,
            is_read1 = record.is_read1,
            cigar = record.cigartuples
        )

    @staticmethod
    def __determine_strand__(record: pysam.AlginmentSegment) -> str:
        u"""
        determine the strand from pysam.AlignmentSegment
        :param record: pysam.AlignmentSegment
        """

        if record.is_read1:
            return "-" if record.is_reverse else "+"
        
        if record.is_read2:
            return "+" if record.is_reverse else "-"

        if not record.is_paired:
            return "-" if record.is_reverse else "+"
        return "*"

    @property
    def exon_parts(self) -> List[int]:
        u"""
        generate exon parts from cigar

        M	BAM_CMATCH	0
        I	BAM_CINS	1
        D	BAM_CDEL	2
        N	BAM_CREF_SKIP	3
        S	BAM_CSOFT_CLIP	4
        H	BAM_CHARD_CLIP	5
        P	BAM_CPAD	6
        =	BAM_CEQUAL	7
        X	BAM_CDIFF	8
        B	BAM_CBACK	9

        :return list of int
        """
        pos = self.start
        pos_list = [pos]

        skipped_code = (1, 2, 4, 5)
        
        for i, j in record.cigartuples:

            if i not in skipped_code:
                pos += j

            if i == 3:
                pos_list.append(pos - j)
                pos_list.append(pos - 1)
        
        pos_list.append(self.end)

        return pos_list



if __name__ == '__main__':
    print(Region.create_from_bed("chr1\t1\t22\t.\t.\t+\tasfd"))
    pass