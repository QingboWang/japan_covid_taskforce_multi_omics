import pandas as pd

#1. Severity stats per batch
df = pd.read_csv("~/Downloads/n1404_parsed_matrix_invnorm.tsv.gz", sep='\t', compression='gzip', index_col=0, header=[0,1])
cov = pd.read_csv("~/Downloads/n1384.protein.combined_covariates_idadded.txt", sep="\t", index_col=0) #moved from shirokane
cov = cov.T #row = samples, col = covariates. For Xs too, row = samples
df = df.T
df.reset_index(level="batch", inplace=True)
cov = cov.join(df.batch, how="left")
st = cov.groupby(["batch", "severity"]).size().unstack().fillna(0).astype(int)
st.columns = ["Asymptomatic", "Mild", "Severe", "Most Severe"]
st.columns.name = "Severity"
frac = (st.T / st.sum(axis=1)).T
err = np.sqrt((frac*(1-frac)).T/st.sum(axis=1)).T
#Plot:
colors = ["tab:gray", "lightsteelblue", "royalblue", "blueviolet"]
frac.plot.bar(stacked=True, color=colors, figsize=(6,3.5)).legend(loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.xticks([0,1,2], ["Batch 1", "Batch 2", "Batch 3"], rotation=30)
plt.xlabel("Protein measurement batch")
plt.ylabel("Sample severity fraction")
for i in range(st.shape[0]):
    plt.text(x=i - 0.49, y=0.8, s="n=")
    for j in range(st.shape[1]):
        plt.text(x=i-0.5, y=0.02+0.22*j, s=st.iloc[i,j], color=colors[j])
plt.ylim([0,1])
plt.tight_layout()
plt.savefig("/Users/qingbowang/Downloads/severity_per_batch.pdf", dpi=500)
plt.savefig("/Users/qingbowang/Downloads/severity_per_batch.png", dpi=500)
plt.clf()


#Check the guys below LOD
ol1filt = pd.read_csv("~/Downloads/Q-00984_NPX_idannot_filtered_to_CT.csv", sep=';')
ol2filt = pd.read_csv("~/Downloads/Q-01085_NPX_idannot_filtered_to_CT.csv", sep=';')
ol3filt = pd.read_csv("~/Downloads/Q-01086_NPX_idannot_filtered_to_CT.csv", sep=';')
#How is LOD defined, per batch or per gene? Check e.g. taking CLEC4C as an example:
ol1filt.LOD.hist(bins=100)
plt.show()#Not consistent of course
ol1filt[ol1filt.Assay=="CLEC4C"].LOD.hist(bins=100)
plt.show()#Very close to constant..??
ol1filt[ol1filt.Assay=="CLEC4C"].LOD.sort_values()
ol1filt[ol1filt.Assay=="CLEC4C"].LOD.value_counts() #Oh OK this is defined per plate

#Just filter to the ones below LOD:
ol1filt["below_lod"] = ol1filt.NPX < ol1filt.LOD
ol1filt["below_lod"].value_counts() #Hmm exists yes...
tb1 = ol1filt.drop_duplicates(subset=["SubjectID","Assay"], keep="first").set_index(["SubjectID","Assay"]).below_lod.unstack()
ol2filt["below_lod"] = ol2filt.NPX < ol2filt.LOD
ol2filt["below_lod"].value_counts() #Hmm exists yes...
tb2 = ol2filt.drop_duplicates(subset=["Variable 2 (Subject)","Assay"], keep="first").set_index(["Variable 2 (Subject)","Assay"]).below_lod.unstack()
ol3filt["below_lod"] = ol3filt.NPX < ol3filt.LOD
ol3filt["below_lod"].value_counts() #Hmm exists yes...
tb3 = ol3filt.drop_duplicates(subset=["Variable 2 (Subject)","Assay"], keep="first").set_index(["Variable 2 (Subject)","Assay"]).below_lod.unstack()
tb2.index.name = tb1.index.name
tb3.index.name = tb1.index.name
tb = pd.concat([tb1, tb2, tb3], axis=0)
tb.to_csv("~/Downloads/olink_flag_below_lod.tsv.gz", sep='\t', compression="gzip")

lodflag = pd.read_csv("~/Downloads/olink_flag_below_lod.tsv.gz", sep='\t', compression="gzip", index_col=0)
#mat = pd.read_csv("~/Downloads/n1384.protein.expression.bed", sep="\t", index_col=[0,1,2,3])
#ipQTL stats:
p_cov_ageon = pd.read_csv("~/Downloads/ip_cis_qtl_top_assoc_reimputed.txt", sep="\t", index_col=0)
e_cov_ageon = pd.read_csv("~/Downloads/ie_cis_qtl_top_assoc_reimputed.txt", sep="\t", index_col=0)
p_cov_ageon = p_cov_ageon.join(e_cov_ageon.pval_gi, how="left", rsuffix="_e")
#gene name
ensgids = pd.read_csv("~/Downloads/gencode.v30.genes.parsed.tsv.gz", sep='\t')
ensgids["ensgid_unv"] = ensgids.ensg_id.str.split("\\.").str[0]
ensgids.index =ensgids.ensgid_unv
p_cov_ageon.index = p_cov_ageon.index.str.split("\\.").str[0]
p_cov_ageon = p_cov_ageon.join(ensgids.gene_name, how="left")
p_cov_ageon.reset_index(inplace=True)
p_cov_ageon.set_index("gene_name", inplace=True)
lodfrac = lodflag.mean(axis=0) #Intersecting samples for bridge-norm are counted twice
lodfrac.name = "perc_below_lod"
p_cov_ageon = p_cov_ageon.join(lodfrac, how="left")
import matplotlib.cm as cm
p_cov_ageon.sort_values(by="perc_below_lod", ascending=True, inplace=True)
plt.figure(figsize=(7,3.5))
plt.scatter(-np.log10(p_cov_ageon.pval_gi_e), -np.log10(p_cov_ageon.pval_gi), c=p_cov_ageon.perc_below_lod*100, cmap='viridis')
plt.plot([0,35], [0,35], linestyle="--", linewidth=1, color="black", zorder=-2)
plt.xlim([0,35])
plt.ylim([0,35/2])
plt.xlabel("Interaction eQTL -log10(p)")
plt.ylabel("Interaction pQTL -log10(p)")
plt.colorbar(label="%Protein measurement below LOD")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Downloads/ie_vs_ipQTL_perlod.png", dpi=500)
plt.savefig("/Users/qingbowang/Downloads/ie_vs_ipQTL_perlod.pdf", dpi=500)
plt.clf()
#For manual annotation:
p_cov_ageon[p_cov_ageon.perc_below_lod>0.1].sort_values(by="pval_gi_e", ascending=True)
p_cov_ageon.sort_values(by="pval_gi_e", ascending=True)