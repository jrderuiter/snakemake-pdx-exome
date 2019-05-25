from os import path

import numpy as np
import pandas as pd


rule sort_targets:
    input:
        config["freebayes"]["targets"]
    output:
        "vcf/split/targets.sorted.bed"
    shell:
        "sort -k1,1 -k2,2n {input[0]} > {output[0]}"


rule merge_targets:
    input:
        "vcf/split/targets.sorted.bed"
    output:
        "vcf/split/targets.merged.bed"
    conda:
        path.join(workflow.basedir, "envs/bedtools.yaml")
    shell:
        "bedtools merge -i {input[0]} > {output[0]}"


rule split_targets:
    input:
        "vcf/split/targets.merged.bed"
    params:
        parts=config["freebayes"]["threads"]
    output:
        ["vcf/split/targets.{}.bed".format(i + 1)
         for i in range(config["freebayes"]["threads"])]
    run:
        targets = pd.read_csv(input[0], sep="\t", header=None)

        if len(targets) < params.parts:
            raise ValueError("Number of freebayes threads must be less than or "
                             "equal to the number of (merged) target regions")

        chunk_size = len(targets) / params.parts
        chunks = (np.arange(len(targets)) / chunk_size).astype(int)

        for i, chunk in targets.groupby(chunks):
            chunk.to_csv(output[i], sep="\t", index=False, header=None)


rule freebayes:
    input:
        samples=expand("bam/deduped/{sample}.bam", sample=get_samples()),
        indices=expand("bam/deduped/{sample}.bam.bai", sample=get_samples()),
        targets="vcf/split/targets.{part}.bed"
    output:
        "vcf/split/calls.{part}.vcf"
    params:
        reference=config["freebayes"]["reference"],
        extra=config["freebayes"]["extra"]
    log:
        "logs/freebayes.log"
    conda:
        path.join(workflow.basedir, "envs/freebayes.yaml")
    shell:
        "freebayes {params.extra} --fasta-reference {params.reference}"
        " --targets {input.targets} {input.samples} > {output[0]} 2> {log}"


# TODO: Sort vcf?

rule vcf_merge:
    input:
        expand("vcf/split/calls.{part}.vcf",
               part=[i + 1 for i in range(config["freebayes"]["threads"])])
    params:
        "--multiallelics -any"
    output:
        "vcf/calls.vcf"
    threads:
        config["bcftools"]["threads"]
    conda:
        path.join(workflow.basedir, "envs/bcftools.yaml")
    shell:
        "bcftools concat --allow-overlaps --remove-duplicates"
        " {params} -o {output[0]} --threads {threads} {input[0]}"


rule vcf_clean:
    input:
        "vcf/calls.vcf"
    params:
        reference="/path/to/reference"
    output:
        "vcf/calls.clean.vcf"
    conda:
        path.join(workflow.basedir, "envs/vcf_utilities.yaml")
    shell:
        # TODO: Fix ambig?
        # TODO: Remove missing alt?
        "cat {input[0]} | "
        "bcftools filter -i 'ALT=\"<*>\" || QUAL > 5' | "
        "vcfallelicprimitives -t DECOMPOSED --keep-geno | "
        "vcffixup - | "
        "vcfstreamsort | "
        "vt normalize -n -r {params.reference} -q - | "
        "vcfuniqalleles | "
        "vt uniq > {output[0]}"

        # Subset alleles for targets (after filter).
        # bcftools view {samples} -a -

# bcftools filter -i 'ALT="<*>" || QUAL > 5' """
# "| {fix_ambig} | bcftools view {samples} -a - | "
# "{py_cl} -x 'bcbio.variation.freebayes.remove_missingalt(x)' | "
# "vcfallelicprimitives -t DECOMPOSED --keep-geno | vcffixup -
# | vcfstreamsort | " "vt normalize -n -r {ref_file} -q -
# | vcfuniqalleles | vt uniq


# rule vcf_soft_filter:


# rule extract_germline:
#     input:
#         "vcf/calls.vcf"
#     output:
#         "vcf/calls.germline.vcf"
#     params:
#         general_filter="",
#         normal_filter="isVariant(GEN[{sample}])",
#         normals=["S1", "S2"],
#         extra=""
#     conda:
#         path.join(workflow.basedir, "envs/snpsift.yaml")
#     script:
#         path.join(workflow.basedir, "scripts/snpsift_extract_germline.py")


# rule annotate_germline:
#     input:
#         calls="vcf/calls.vcf",
#         germline_calls="vcf/calls.germline.vcf"
#     output:
#         "vcf/calls.somatic.vcf"
#     conda:
#         path.join(workflow.basedir, "envs/snpsift.yaml")
#     shell:
#         "snpsift annotate {params} -exists Germline"
#         " {input.germline_calls} {input.calls} > {output[0]}"


# rule annotate_database:
#     input:
#         calls=
#         database=
#     output:

#     shell:
#         "snpsift annotate {params} -id"
#         " {input.germline_calls} {input.calls} > {output[0]}"

# ftp://ftp-mouse.sanger.ac.uk//REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
# ftp://ftp-mouse.sanger.ac.uk//REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz
