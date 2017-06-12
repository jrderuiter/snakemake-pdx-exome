
rule disambiguate:
    input:
        human='merged/{sample}.human.bam'),
        mouse='merged/{sample}.mouse.bam')
    output:
        temp('disambiguate/{sample}.ambiguousSpeciesA.bam')),
        temp('disambiguate/{sample}.ambiguousSpeciesB.bam')),
        temp('disambiguate/{sample}.disambiguatedSpeciesA.bam')),
        temp('disambiguate/{sample}.disambiguatedSpeciesB.bam'))
    params:
        output_dir='disambiguate',
        algorithm='bwa',
        prefix='{sample}'
    shell:
        'ngs_disambiguate -o {params.output_dir} -s {params.prefix}'
        ' -a {params.algorithm} {input.human} {input.mouse}'


rule disambiguate_sort:
    input:
        'disambiguate/{sample}.disambiguatedSpeciesA.bam')
    output:
        'filtered/{sample}.bam'
    params:
        options=' '.join(config['picard_sort_sam']['options'] or [])
    log:
        'logs/disambiguate_sort/{sample}.log')
    shell:
        'picard SortSam {params.options} INPUT={input}'
        ' OUTPUT={output} SORT_ORDER=coordinate'
        ' VALIDATION_STRINGENCY=LENIENT 2> {log}'
