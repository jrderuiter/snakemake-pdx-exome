from snakemake.shell import shell


# Build filter expression.
normal_filters = [snakemake.params.normal_filter.format(sample)
                  for sample in snakemake.params.normals]

filters = [snakemake.params.general_filter] + normal_filters

filter_expr = " & ".join([f for f in filters if f])

# Run command.
shell("snpsift filter"
      " {snakemake.params.extra}"
      " {filter_expr}"
      " {snakemake.input[0]}"
      " > snakemake.output[0]")
