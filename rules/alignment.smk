from os import path


rule bwa_graft:
    input:
        ['fastq/trimmed/{sample}.{lane}.R1.fastq.gz',
         'fastq/trimmed/{sample}.{lane}.R2.fastq.gz'],
    output:
        temp('bam/aligned/{sample}.{lane}.graft.bam')
    params:
        index=config['bwa_graft']['index'],
        extra=config['bwa_graft']['extra'],
        sort='picard',
        sort_order='queryname',
        sort_extra=config['bwa_graft']['sort_extra']
    threads:
        config['bwa_graft']['threads']
    log:
        'logs/bwa/{sample}.{lane}.graft.log'
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/bwa/mem')


rule bwa_host:
    input:
        ['fastq/trimmed/{sample}.{lane}.R1.fastq.gz',
         'fastq/trimmed/{sample}.{lane}.R2.fastq.gz']
    output:
        temp('bam/aligned/{sample}.{lane}.host.bam')
    params:
        index=config['bwa_host']['index'],
        extra=config['bwa_host']['extra'],
        sort='picard',
        sort_order='queryname',
        sort_extra=config['bwa_host']['sort_extra']
    threads:
        config['bwa_host']['threads']
    log:
        'logs/bwa/{sample}.{lane}.host.log'
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/bwa/mem')


def merge_inputs(wildcards):
    lanes = get_sample_lanes(wildcards.sample)

    file_paths = ['bam/aligned/{}.{}.{}.bam'.format(
                    wildcards.sample, lane, wildcards.organism)
                  for lane in lanes]

    return file_paths


rule picard_merge_bam:
    input:
        merge_inputs
    output:
        'bam/merged/{sample}.{organism}.bam'
    params:
        config['picard_merge_bam']['extra']
    log:
        'logs/picard_merge_bam/{sample}.{organism}.log'
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/picard/mergesamfiles')


rule disambiguate:
    input:
        a='bam/merged/{sample}.graft.bam',
        b='bam/merged/{sample}.host.bam'
    output:
        a_ambiguous=temp('bam/disambiguate/{sample}.graft.ambiguous.bam'),
        b_ambiguous=temp('bam/disambiguate/{sample}.host.ambiguous.bam'),
        a_disambiguated=temp('bam/disambiguate/{sample}.graft.bam'),
        b_disambiguated=temp('bam/disambiguate/{sample}.host.bam'),
        summary='bam/disambiguate/{sample}.summary.txt'
    params:
        algorithm='bwa',
        prefix='{sample}',
        extra=config['disambiguate']['extra']
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/ngs-disambiguate')


rule picard_sort_bam:
    input:
        'bam/disambiguate/{sample}.graft.bam'
    output:
        'bam/sorted/{sample}.bam'
    params:
        sort_order='coordinate',
        extra=config['picard_sort_bam']['extra']
    log:
        'logs/picard_sort_bam/{sample}.log'
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/picard/sortsam')


rule picard_mark_duplicates:
    input:
        'bam/sorted/{sample}.bam'
    output:
        bam='bam/deduped/{sample}.bam',
        metrics='bam/deduped/{sample}.metrics'
    params:
        config['picard_mark_duplicates']['extra']
    log:
        'logs/picard_mark_duplicates/{sample}.log'
    wrapper:
        '0.15.4/bio/picard/markduplicates'


rule samtools_index:
    input:
        'bam/deduped/{sample}.bam'
    output:
        'bam/deduped/{sample}.bam.bai'
    wrapper:
        "0.15.4/bio/samtools/index"
