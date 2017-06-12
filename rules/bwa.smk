
rule bwa_align_human:
    input:
        fastq1='trimmed/{sample}.{lane}.R1.fastq.gz',
        fastq2='trimmed/{sample}.{lane}.R1.fastq.gz'
    output:
        temp('aligned/{sample}.{lane}.human.bam')
    params:
        options=' '.join(config['bwa_align_human']['options'] or []),
        index=config['bwa_align_human']['index'],
        rg=lambda wildcards: readgroup_str(wildcards),
        sort_options=' '.join(config['picard_sort_sam']['options'] or [])
    threads:
        config['bwa_align_human']['threads']
    log:
        bwa='logs/bwa/{sample}.{lane}.human.align.log',
        sort='logs/bwa/{sample}.{lane}.human.sort.log'
    shell:
        'bwa mem -t {threads} {params.options} {params.index} -R "{params.rg}"'
        ' {input.fastq1} {input.fastq2} 2> {log.bwa} | '
        'picard SortSam {params.sort_options} INPUT=/dev/stdin'
        ' OUTPUT={output} SORT_ORDER=queryname 2> {log.sort}'

rule bwa_align_mouse:
    input:
        fastq1='trimmed/{sample}.{lane}.R1.fastq.gz',
        fastq2='trimmed/{sample}.{lane}.R1.fastq.gz'
    output:
        temp('aligned/{sample}.{lane}.human.bam')
    params:
        options=' '.join(config['bwa_align_mouse']['options'] or []),
        index=config['bwa_align_mouse']['index'],
        rg=lambda wildcards: readgroup_str(wildcards),
        sort_options=' '.join(config['picard_sort_sam']['options'] or [])
    threads:
        config['bwa_align_mouse']['threads']
    log:
        bwa='logs/bwa/{sample}.{lane}.mouse.align.log',
        sort='logs/bwa/{sample}.{lane}.mouse.sort.log'
    shell:
        'bwa mem -t {threads} {params.options} {params.index} -R "{params.rg}"'
        ' {input.fastq1} {input.fastq2} 2> {log.bwa} | '
        'picard SortSam {params.sort_options} INPUT=/dev/stdin'
        ' OUTPUT={output} SORT_ORDER=queryname 2> {log.sort}'
