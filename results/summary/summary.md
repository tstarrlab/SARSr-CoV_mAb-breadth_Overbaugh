# Summary

Analysis run by [Snakefile](../../Snakefile)
using [this config file](../../config.yaml).
See the [README in the top directory](../../README.md)
for details.

Here is the DAG of the computational workflow:
![dag.svg](dag.svg)

Here is the Markdown output of each Jupyter notebook in the
workflow:



1. Get prior sarbecovirus homolog wildtypes expression measurements from [this repository](https://github.com/jbloomlab/SARSr-CoV_homolog_survey).

2. [Count variants by barcode](count_variants.md).
   Creates a [variant counts file](../counts/variant_counts.csv)
   giving counts of each barcoded variant in each condition.

3. [Fit EC50 to mAb binding curves](compute_EC50.md).
   Creates a [table](../bc_mAb_EC50/bc_mAb_EC50.csv)
   giving the EC50 phenotype of each barcoded variant in each condition.

4. Collapse internal replicate barcodes of each variant to final variant phenotypes for the wildtype sarbecovirus homologs pool. Analysis [here](collapse_barcodes_lib47_SARSr-wts.md) and final output file [here](../final_variant_scores/final_variant_scores_lib47.csv).