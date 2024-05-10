## Analysis
Codes relevant to the main analysis results are in this `Analysis` directory.
- `ukb_replication.py` was used for replication of our pQTL results in light with UKB-PPP data (Fig. 1)
- `pqtl_analysis_and_plots.py` was used for pQTL analysis (Fig. 1, 2)
- `specificqtls_analysis.py` was used for analysis of mRNA or protein specific QTLs (Fig. 3)
- `classify_proteins.py` was used to analyse proteins in light of blood mRNA expression levels and other classifications (Fig. 3, 6)
- `coloc_analysis.py` was used for colocalization analysis of e/pQTLs (Fig. 3, 4, 7)
- `locuszoom_plots.py` was used to generate locus zoom figures (Fig. 3, 4, 7)
- `trans_qtl_analysis.py` was used for trans-e/pQTL analysis (Fig. 5)
- `eqtl_analysis_and_plots.py` was used for eQTL analysis (Ext. Fig. 2)
- `mpra_plots.py` was used for analysis of MPRA results (Ext. Fig. 2)


## Supple
Codes relevant to the supplementary data are in this `Supple` directory.
- `susie_coloc.py` was used for colocalization analysis using the SuSiE-coloc framework
- `sqtl_call.sh` and `sqtl_analysis.py` were used for sQTL analysis
- `rare_var_inspection.py` was used to evaluate the rare variants' functional enrichment and other properties
- `protein_lod_related.py` was used to investigate the effect of Limitation of Detection (LoD) of protein measurements
- `network_analysis.py` was used to investigate the mRNA and protein expression modules
- `MR_analysis.py` and `MR_analysis.R` was used for testing the mRNA/protein specific regulation using Mendelian Randomization tools
- `eqtlgen_replication.py` was used for comparison with the eQTLgen data
- `covid_hgi_intersection.py` was used to test the association signal overlap with COVID-19 HGI data

## qc
Codes relevant to the quality control / data generation steps are in this `qc` directory.
- `rna_qc.py` was used for quality control of mRNA expression matrix
- `protein_qc_plots.py` was used for quality control of protein expression matrix
- `prep_protein_expression_matrix.py` and `prep_protein_expression_matrix.R` was used to prepare the protein expression matrix (e.g. normalization within and across batch)
- `create_ps_and_qs_matrix.py` was used to create mRNA-adjusted protein expression matrix and protein-adjusted mRNA expression matrix



## fine_mapping
Codes relevant to the fine-mapping pipeline are in this `fine_mapping` directory.
- `create_zs.py` and `parse_sumstats.py` were used to parse the summary statistics from fastQTL into the format compatible for fine-mapping algorithms
- `create_ld_mat_for_a_gene.py` was used to create in-sample LD matrix that was feed into the fine-mapping algorithms
- `finemap_a_gene.sh` was used to fine-map eQTLs for a gene (and es/ps QTLs as well)
- `finemap_a_gene_pqtl.sh` was used to fine-map pQTLs for a protein
- `finemap_sqtls.sh` was used to fine-map sQTLs for a gene
- `run_susie_for_a_gene.R` is the code to run SuSiE fine-mapping
- `collect_fm_results.sh` was used to parse the fine-mapping outputs.

