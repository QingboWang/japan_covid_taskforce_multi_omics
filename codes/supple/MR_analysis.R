#as in https://mrcieu.github.io/TwoSampleMR/articles/introduction.html
#library(remotes)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#First, simply eQTL vs pQTL sumstats:
i=1
for (chr in seq(1,22)){
  e = read.table(paste0("~/Downloads/n998_eqtl_sumstats/eqtl_n998.chr", chr, ".allpairs.txt.gz"),sep='\t', header=TRUE)
  p = read.table(paste0("~/Downloads/n998_pqtl_sumstats/pqtl_n998.chr", chr, ".allpairs.txt.gz"),sep='\t', header=TRUE)
  refminor = read.table(paste0("~/Downloads/n998_refminor/n998_ref_is_minor_chr", chr, ".tsv.gz"),sep='\t', header=TRUE)
  refminor$ref_is_minor = TRUE
  e = e[e$ma_count>2,]
  p = p[p$ma_count>2,]

  #annotate alt_af etc
  library(dplyr)
  refminor$variant_id = refminor$ID
  e = left_join(e, refminor[,c("variant_id", "ref_is_minor")], by = "variant_id")
  e$ref_is_minor = ifelse(is.na(e$ref_is_minor), FALSE, e$ref_is_minor)
  table(e$ref_is_minor)#sanity check: Looks good (なぜ今回は半分程度にならんのだ？??? )
  p = left_join(p, refminor[,c("variant_id", "ref_is_minor")], by = "variant_id")
  p$ref_is_minor = ifelse(is.na(p$ref_is_minor), FALSE, p$ref_is_minor)

  #annotate ref, alt and alt_af
  e$variant_id_hg38 = sapply(strsplit(e$variant_id, "_"), function(x) x[1])
  e$ref = sapply(strsplit(e$variant_id_hg38, ":"), function(x) x[3])
  e$alt = sapply(strsplit(e$variant_id_hg38, ":"), function(x) x[4])
  p$variant_id_hg38 = sapply(strsplit(p$variant_id, "_"), function(x) x[1])
  p$ref = sapply(strsplit(p$variant_id_hg38, ":"), function(x) x[3])
  p$alt = sapply(strsplit(p$variant_id_hg38, ":"), function(x) x[4])

  e$alt_af = ifelse(e$ref_is_minor, 1-e$maf, e$maf)
  p$alt_af = ifelse(p$ref_is_minor, 1-p$maf, p$maf)

  e$qtl_class = "eqtl"
  p$qtl_class = "pqtl"
  e$qtl_id = "eqtl_jctf"
  p$qtl_id = "pqtl_jctf"

  #make sure that the order is the same
  #table(paste0(e$gene_id, e$variant_id)==paste0(p$gene_id, p$variant_id)) #perfect.

  #iterate through genes
  genes = unique(e$gene_id)
  for (gn in genes){
    e1 = e[e$gene_id==gn,]
    p1 = p[p$gene_id==gn,]

    # List available GWASs
    #ao <- available_outcomes()
    # Get instruments
    #exposure_dat <- extract_instruments("ieu-a-2") #These lines are when from reference

    exposure_dat = e1[c("variant_id", "slope", "slope_se", "alt", "ref", "alt_af", "qtl_class", "qtl_id")]
    cols_exposure = c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure",
      "other_allele.exposure", "eaf.exposure", "exposure", "id.exposure")
    colnames(exposure_dat) = cols_exposure

    outcome_dat = p1[c("variant_id", "slope", "slope_se", "alt", "ref", "alt_af", "qtl_class", "qtl_id")]
    cols_outcome = c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome",
      "other_allele.outcome", "eaf.outcome", "outcome", "id.outcome")
    colnames(outcome_dat) = cols_outcome

    # Harmonise the exposure and outcome data (which is already harmonized by definition, actually)
    dat = harmonise_data(exposure_dat, outcome_dat)
    # Perform MR
    e_to_p = mr(dat)
    e_to_p$gene_id = gn

    #Also do the opposite
    exposure_dat2 = data.frame(outcome_dat)
    outcome_dat2 <- data.frame(exposure_dat)
    colnames(outcome_dat2) = cols_outcome
    colnames(exposure_dat2) = cols_exposure
    dat2 = harmonise_data(exposure_dat2, outcome_dat2)
    p_to_e = mr(dat2)
    p_to_e$gene_id = gn

    write.table(e_to_p, paste0("~/Downloads/mr_results/", gn, "_e_to_p.tsv"), sep="\t", quote=FALSE)
    write.table(p_to_e, paste0("~/Downloads/mr_results/", gn, "_p_to_e.tsv"), sep="\t", quote=FALSE)
    print (paste0("Done ",i, ", ", Sys.time()))
    i = i+1
  }
}





