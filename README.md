# scATS

ATS idenfication and quantification tools

>This code still under development

## Installation

```bash
git clone git@github.com:ygidtu/scATS.git

cd scATS

# install scATS as command line tool
python3 setup.py install

# Or just run source code
python3 main.py
```

## Usage

```bash
Usage: main.py [OPTIONS] COMMAND [ARGS]...

  Welcome

  This function is used to test the function of sashimi plotting

  Created by ygidtu@gmail.com at 2018.12.19 :return:

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  ats      Inference :param debug: enable debug mode
  count    Postprocess: count ats and calculate psi
```

### ats

Identify ATSs based on aligned BAM files.

```bash
➜  afe git:(master) ✗ python main.py ats --help
Usage: main.py ats [OPTIONS] BAMS...

  Inference

  :param debug: enable debug mode

Options:
  -g, --gtf PATH                 The path to genome annotation file in GTF format.   [required]
  -o, --output PATH              The path to output file.   [required]
  -u, --utr-length INTEGER       The length of UTR.
  --n-max-ats INTEGER RANGE      The maximum number of ATSs in same UTR.
  --n-min-ats INTEGER RANGE      The minimum number of ATSs in same UTR.
  --min-ws FLOAT                 The minimum weight of ATSs.
  --min-reads INTEGER            The minimum number of reads in UTR.
  --max-unif-ws FLOAT            The maximum weight of uniform component.
  --max-beta INTEGER             The maximum std for ATSs.
  --fixed-inference              Inference with fixed parameters.
  -d, --debug                    Enable debug mode to get more debugging
                                 information, Never used this while
                                 running.
  -p, --processes INTEGER RANGE  How many cpu to use.
  --remove-duplicate-umi         Only kept reads with different UMIs for
                                 ATS inference.
  -h, --help                     Show this message and exit.
```

### count

Quantification of ATSs.

```bash
➜  afe git:(master) ✗ python main.py count --help
Usage: main.py count [OPTIONS]

  Postprocess: count ats and calculate psi

Options:
  -i, --ats PATH                 The path to inferred ats sites.
  -b, --bam PATH                 The file contains path to bams.
  --delimiter TEXT               The delimiter of input bam list
  -o, --output PATH              The path to output utr file, bed format.
  -p, --processes INTEGER RANGE  How many cpu to use.
  -c, --compress                 Wheter to save in gzip format.
  --bulk                         Wheter the input bam is Nanopore or
                                 PacBio.
  -h, --help                     Show this message and exit.
```

## Example

Please prepare the genome annotation file in GTF format and corresponding genome sequance file in fasta format first.

```bash
cd simulation

# generate simulation data
python simulation.py /path/to/gtf /path/to/fasta --n_jobs 2 --total 10000

# run ats model
python ../main.py ats -g ./tss.gtf -o ./inferred_sites.txt -p 4 ./simu.bam
```

### The output file details

```bash
utr     gene_name       transcript_name number_of_reads inferred_sites   alpha_arr       beta_arr        ws      L
1:1168755-1169255:+     MIR429  MIR429-201      27      1168924,1169077,1169096,1169101 169,322,341,346 5,5,10,5        0.07156744317872403,0.8035226895855316,0.03289719913137757,0.08391005080588618,0.008102617298480465     500
1:1891221-1891721:+     AL109917.1      AL109917.1-201  10      1891223,1891474,1891539 2,253,318       5,5,10  0.3979824461238993,0.2978253162745582,0.2929895754731792,0.01120266212836311    500
```

- utr: the genomic location of UTR.
- gene_name: the name of host gene of this UTR.
- transcript_name: the name of host transcripts of this UTR.
- number_of_reads: the number of reads used to infer ATS sites in corresponding UTR.
- inferred_sites: the inferred ATS sites, multiple sites were seperate by comma.
- alpha_arr: the alpha of guassion distribution.
- beta_arr: the beta of guassion distribution.
- ws: weights of each ATS sites, and weights of reads not belong to each ATS sites.
- L: the length of UTR, which may larger than 500 (default UTR length). Due to overlapped UTR merging.
