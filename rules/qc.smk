from os import path


rule multiqc:
    input:
        directory='.',
        fastqc=expand('qc/fastqc/{sample_lane}.{pair}_fastqc.html',
                      sample_lane=get_samples_with_lane(), pair=['R1', 'R2']),
        samtools_stats=expand('qc/samtools_stats/{sample}.txt',
                              sample=get_samples())
    output:
        'qc/multiqc_report.html'
    params:
        output_dir='qc',
        extra=''
    shell:
        'multiqc {params.extra} --force -o'
        ' {params.output_dir} {input.directory}'


rule fastqc:
    input:
        'fastq/trimmed/{sample}.{lane}.{pair}.fastq.gz'
    output:
        html='qc/fastqc/{sample}.{lane}.{pair}_fastqc.html',
        zip='qc/fastqc/{sample}.{lane}.{pair}_fastqc.zip'
    params:
        config['fastqc']['extra']
    wrapper:
        'file://' + path.join(workflow.basedir, 'wrappers/fastqc')


rule samtools_stats:
    input:
        'bam/deduped/{sample}.bam'
    output:
        'qc/samtools_stats/{sample}.txt'
    shell:
        'samtools stats {input} > {output}'
