import pandas as pd

configfile: 'config.yaml'

################################################################################
# Globals                                                                      #
################################################################################

samples = pd.read_csv('samples.tsv', sep='\t')


################################################################################
# Functions                                                                     #
################################################################################

def get_samples():
    return list(samples['sample'].unique())

def get_sample_lanes(sample):
    subset = samples.loc[samples['sample'] == sample]
    return list(subset['lane'].unique())


################################################################################
# Rules                                                                        #
################################################################################

rule all:
    input:
        #expand('filtered/{sample}.bam', sample=get_samples())
        'fastq/1957_6_T250.L001.R1.fastq.gz'

include: "rules/fastq.smk"
include: "rules/cutadapt.smk"
include: "rules/fastqc.smk"
include: "rules/bwa.smk"
include: "rules/disambiguate.smk"
