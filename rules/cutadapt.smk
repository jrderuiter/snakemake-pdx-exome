
rule cutadapt:
    input:
        'fastq/{sample}.{pair}.fastq.gz'
    output:
        fastq=temp('trimmed/{sample}.{pair}.fastq.gz'),
        qc='trimmed/{sample}.{pair}.qc.txt'
    params:
        options=config['cutadapt']['options']
    log:
        path.join('logs/cutadapt/{sample}.{pair}.log'),
    shell:
        'cutadapt {params.options} -o {output.fastq} {input}'
        ' > {output.qc} 2> {log}'
