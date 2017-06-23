from os import path


rule bwa_graft:
    input:
        ["fastq/trimmed/{sample}.{lane}.R1.fastq.gz",
         "fastq/trimmed/{sample}.{lane}.R2.fastq.gz"],
    output:
        temp("bam/aligned/{sample}.{lane}.graft.bam")
    params:
        index=config["bwa"]["index_graft"],
        extra=config["bwa"]["extra"],
        sort="samtools",
        sort_order="queryname",
        sort_extra=config["bwa"]["sort_extra"]
    threads:
        config["bwa"]["threads"]
    log:
        "logs/bwa/{sample}.{lane}.graft.log"
    wrapper:
        "master/bio/bwa/mem"


rule bwa_host:
    input:
        ["fastq/trimmed/{sample}.{lane}.R1.fastq.gz",
         "fastq/trimmed/{sample}.{lane}.R2.fastq.gz"]
    output:
        temp("bam/aligned/{sample}.{lane}.host.bam")
    params:
        index=config["bwa"]["index_host"],
        extra=config["bwa"]["extra"],
        sort="samtools",
        sort_order="queryname",
        sort_extra=config["bwa"]["sort_extra"]
    threads:
        config["bwa"]["threads"]
    log:
        "logs/bwa/{sample}.{lane}.host.log"
    wrapper:
        "master/bio/bwa/mem"


rule disambiguate:
    input:
        a="bam/aligned/{sample}.{lane}.graft.bam",
        b="bam/aligned/{sample}.{lane}.host.bam"
    output:
        a_ambiguous=temp("bam/disambiguate/{sample}.{lane}.graft.ambiguous.bam"),
        b_ambiguous=temp("bam/disambiguate/{sample}.{lane}.host.ambiguous.bam"),
        a_disambiguated=temp("bam/disambiguate/{sample}.{lane}.graft.bam"),
        b_disambiguated=temp("bam/disambiguate/{sample}.{lane}.host.bam"),
        summary="qc/disambiguate/{sample}.{lane}.summary.txt"
    params:
        algorithm="bwa",
        prefix="{sample}.{lane}",
        extra=config["disambiguate"]["extra"]
    wrapper:
        "file://" + path.join(workflow.basedir, "wrappers/ngs-disambiguate")


rule sambamba_sort:
    input:
        "bam/disambiguate/{sample}.{lane}.graft.bam"
    output:
        "bam/sorted/{sample}.{lane}.bam"
    params:
        "--tmpdir=tmp"
    threads: 10
    wrapper:
        "master/bio/sambamba/sort"


def merge_inputs(wildcards):
    lanes = get_sample_lanes(wildcards.sample)

    file_paths = ["bam/sorted/{}.{}.bam".format(wildcards.sample, lane)
                  for lane in lanes]

    return file_paths


rule picard_merge_bam:
    input:
        merge_inputs
    output:
        "bam/merged/{sample}.bam"
    params:
        config["picard_merge_bam"]["extra"]
    log:
        "logs/picard_merge_bam/{sample}.log"
    wrapper:
        "master/bio/picard/mergesamfiles"


rule picard_mark_duplicates:
    input:
        "bam/merged/{sample}.bam"
    output:
        bam="bam/deduped/{sample}.bam",
        metrics="qc/picard_mark_duplicates/{sample}.metrics"
    params:
        config["picard_mark_duplicates"]["extra"]
    log:
        "logs/picard_mark_duplicates/{sample}.log"
    wrapper:
        "0.15.4/bio/picard/markduplicates"


rule samtools_index:
    input:
        "bam/deduped/{sample}.bam"
    output:
        "bam/deduped/{sample}.bam.bai"
    wrapper:
        "0.15.4/bio/samtools/index"
