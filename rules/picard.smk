
def merge_inputs(wildcards):
    lanes = get_sample_lanes(wildcards.sample)

    file_paths = ['aligned/{}.{}.{}.bam'.format(wildcards.sample,
                                                lane,
                                                wildcards.organism)
                  for lane in lanes]

    return file_paths


rule picard_merge_bam:
    input:
        merge_inputs
    output:
        'merged/{sample}.{organism}.bam'
    params:
        options=' '.join(config['picard_merge_bam']['options'] or [])
    log:
        'logs/picard_merge_bam/{sample}.{organism}.log')
    run:
        shell('picard MergeSamFiles {params.options} SORT_ORDER=queryname ' +
              ' '.join('INPUT={}'.format(in_) for in_ in input) +
              ' OUTPUT={output} 2> {log}')
