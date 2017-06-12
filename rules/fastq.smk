
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

def fetch_url(wildcards):
    samples = samples.set_index(['sample', 'lane'])
    pair = 'fastq1' if pair == 'R1' else 'fastq2'
    return samples.loc[(wildcards.sample, wildcards.lane)][pair]

rule fetch:
    input:
        lambda wc: HTTP.remote(fetch_url(wc))
    output:
        'fastq/{sample}.{lane}.{pair}.fastq.gz'
    shell:
        'mv {input} {output}'
