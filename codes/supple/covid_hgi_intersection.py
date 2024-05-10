#Simply look at COVID-19 HGI intersection

import pandas as pd

hgi = pd.read_excel("~/Downloads/41586_2023_6355_MOESM4_ESM.xlsx", sheet_name=1, skiprows=2)

#our JCTF ie and ip QTL results (lead variant):
e = pd.read_csv("~/Downloads/ie_cis_qtl_top_assoc_reimputed.txt", sep="\t", index_col=0)
p = pd.read_csv("~/Downloads/ip_cis_qtl_top_assoc_reimputed.txt", sep="\t", index_col=0)
#map to their gene names:
ensgids = pd.read_csv("~/Downloads/gencode.v30.genes.parsed.tsv", sep='\t')
ensgids["ensg_unversioned"] = ensgids.ensg_id.str.split("\\.").str[0]
ensgids.set_index("gene_name", inplace=True)

coding = hgi.Coding[~hgi.Coding.isna()].str.split(',').explode().reset_index(drop=True).unique()
distance = hgi.Distance[~hgi.Distance.isna()].str.split(',').explode().reset_index(drop=True).unique()
ld = hgi.LD[~hgi.LD.isna()].str.split(',').explode().reset_index(drop=True).unique()
eg = hgi.eGenes[~hgi.eGenes.isna()].str.split(',').explode().reset_index(drop=True).unique()
covid_genes_tier1 = coding
covid_genes_tier2 = np.setdiff1d(np.union1d(distance, np.union1d(ld, eg)), coding)
covid_genes_tier1_ensg = ensgids.loc[ensgids.index.intersection(covid_genes_tier1),"ensg_unversioned"]
covid_genes_tier2_ensg = ensgids.loc[ensgids.index.intersection(covid_genes_tier2),"ensg_unversioned"]

e.index = e.index.str.split("\\.").str[0]
p.index = p.index.str.split("\\.").str[0]

e_risk1 = e.loc[e.index.intersection(covid_genes_tier1_ensg),:]
e_risk2 = e.loc[e.index.intersection(covid_genes_tier2_ensg),:]
e_norisk = e.loc[e.index.difference(covid_genes_tier1_ensg).difference(covid_genes_tier2_ensg),:]
p_risk1 = p.loc[p.index.intersection(covid_genes_tier1_ensg),:]
p_risk2 = p.loc[p.index.intersection(covid_genes_tier2_ensg),:]
p_norisk = p.loc[p.index.difference(covid_genes_tier1_ensg).difference(covid_genes_tier2_ensg),:]

plt.figure(figsize=(5,3.5))
plt.scatter(-np.log10(e_norisk.pval_gi), -np.log10(p_norisk.pval_gi),
            color="tab:blue", alpha=0.75, edgecolor="black", label="No")
plt.scatter(-np.log10(e_risk2.pval_gi), -np.log10(p_risk2.pval_gi),
            color="tab:olive", alpha=0.75, edgecolor="black", label="Non-coding")
plt.scatter(-np.log10(e_risk1.pval_gi), -np.log10(p_risk1.pval_gi),
            color="tab:orange", alpha=0.75, edgecolor="black", label="Coding")
plt.xlabel("COVID-19 interaction eQTL -log10(p)")
plt.ylabel("COVID-19 interaction\npQTL -log10(p)")
plt.legend(title="COVID-19 HGI risk gene", loc="upper center")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Downloads/ieQTL_hgi.png", dpi=500)
plt.savefig("/Users/qingbowang/Downloads/ieQTL_hgi.pdf", dpi=500)
plt.clf()
#And manual annotate probably?
ensgids.reset_index(inplace=True)
ensgids.set_index("ensg_unversioned", inplace=True)
e_risk1 = e_risk1.join(ensgids["gene_name"], how="left")
e_risk1[["pval_gi","gene_name"]].sort_values(by="pval_gi", ascending=True)
p_risk1 = p_risk1.join(ensgids["gene_name"], how="left")
p_risk1[["pval_gi","gene_name"]].sort_values(by="pval_gi", ascending=True)
e_risk2 = e_risk2.join(ensgids["gene_name"], how="left")
e_risk2[["pval_gi","gene_name"]].sort_values(by="pval_gi", ascending=True)
p_risk2 = p_risk2.join(ensgids["gene_name"], how="left")
p_risk2[["pval_gi","gene_name"]].sort_values(by="pval_gi", ascending=True)


