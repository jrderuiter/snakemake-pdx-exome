# Snakemake workflow: pdx-exome

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/pdx-exome.svg?branch=master)](https://travis-ci.org/snakemake-workflows/pdx-exome)

This is a Snakemake workflow for generating variant calls from xenograft exome
sequencing data (or similar targeted DNA-sequencing data from xenograft
samples). The workflow is designed to handle paired-end (and optionally
multi-lane) sequencing data.

The workflow essentially performs the following steps:

* The input reads are trimmed to remove adapters and/or poor quality base calls
  using cutadapt.
* The trimmed reads are aligned to the reference genome using bwa mem.
* The alignments are sorted by queryname using picard SortSam.
* Bam files from multiple lanes are merged using picard MergeSamFiles.
* NGS-disambiguate is used to separate separate graft/host reads.
* The filtered graft bam file is coordinate-sorted using picard SortSam.
* Picard MarkDuplicates is used to remove optical/PCR duplicates.
* Variant calls are generated using freebayes.

QC statistics are generated using fastqc and samtools stats, and are summarized
using multiqc.

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/pdx-exome/releases).
If you intend to modify and further develop this workflow, fork this reposity. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.

## Authors

* Julian de Ruiter (@jrderuiter)
