# Snakemake workflow: pdx-exome

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This is a Snakemake workflow for generating variant calls from xenograft exome
sequencing data (or similar targeted DNA-sequencing data from xenograft
samples). The workflow is designed to handle paired-end (and optionally
multi-lane) sequencing data.

The workflow essentially performs the following steps:

* The input reads are trimmed to remove adapters and/or poor quality base calls
  using cutadapt.
* The trimmed reads are aligned to the reference genome using bwa mem.
  The alignments are sorted by queryname using samtools.
* NGS-disambiguate is used to separate separate graft/host reads.
* The filtered graft bam file is coordinate-sorted using sambamba.
* Bam files from multiple lanes are merged using picard MergeSamFiles.
* Picard MarkDuplicates is used to remove optical/PCR duplicates.
* Variant calls are generated using freebayes.

QC statistics are generated using fastqc and samtools stats, and are summarized
using multiqc.

**Note that this workflow is still under active development.**

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the
[latest release](https://github.com/snakemake-workflows/pdx-exome/releases).
If you intend to modify and further develop this workflow, fork this
repository. Please consider providing any generally applicable modifications
via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits
to the authors by citing the URL of this repository and, if available, its
DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs by editing the config file
`config.yaml` and the samples file `samples.tsv`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run using

    snakemake -n

Execute the workflow locally using

    snakemake --cores $N

using `$N` cores or run it in a cluster environment using

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for
further details.

## Authors

* Julian de Ruiter (@jrderuiter)
