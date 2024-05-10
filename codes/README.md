## Analysis
Codes relevant to the main analysis results are in this `Analysis` directory.

## Supple
Codes relevant to the supplementary data are in this `Supple` directory.

## qc
Codes relevant to the quality control / data generation steps are in this `qc` directory.

## fine_mapping
Codes relevant to the fine-mapping pipeline are in this `fine_mapping` directory.
- `create_zs.py` and `parse_sumstats.py` were used to parse the summary statistics from fastQTL into the format compatible for fine-mapping algorithms
- `create_ld_mat_for_a_gene.py` was used to create in-sample LD matrix that was feed into the fine-mapping algorithms
- `finemap_a_gene.sh` was used to fine-map eQTLs for a gene (and es/ps QTLs as well)
- `finemap_a_gene_pqtl.sh` was used to fine-map pQTLs for a protein
- `finemap_sqtls.sh` was used to fine-map sQTLs for a gene
- `run_susie_for_a_gene.R` is the code to run SuSiE fine-mapping
- `collect_fm_results.sh` was used to parse the fine-mapping outputs.

