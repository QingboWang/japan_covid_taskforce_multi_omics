#collect MR results that was run in R:
import glob
genes = pd.read_csv("~/Downloads/mrna_and_protein_n998.tsv.gz", sep='\t', index_col=0).columns
genes = genes[:int(len(genes)/2)] #removing the _protein guys
e_to_p = []
p_to_e = []
for gn in genes:
    try:
        e_to_p.append(pd.read_csv("/Users/qingbowang/Downloads/mr_results/{0}_e_to_p.tsv".format(gn), sep="\t"))
        p_to_e.append(pd.read_csv("/Users/qingbowang/Downloads/mr_results/{0}_p_to_e.tsv".format(gn), sep="\t"))
    except:
        pass #just that it is not done yet
e_to_p = pd.concat(e_to_p)
p_to_e = pd.concat(p_to_e)
#e.g. take MR Egger example:
d1 = e_to_p[e_to_p.method=="MR Egger"]
d2 = p_to_e[p_to_e.method=="MR Egger"]
repetitions = 100
shuffled_beta_diff = []
for i in range(repetitions):
    result = d1.b.sample(frac=1, random_state=i) - d2.b.sample(frac=1, random_state=i**2)
    shuffled_beta_diff.append(result)
shuffled_beta_diff = pd.concat(shuffled_beta_diff)
real_beta_diff = d1.b - d2.b
real_beta_diff.index = d1.gene_id
def get_pval_betadiff(b):
    return (sum(abs(b)<abs(shuffled_beta_diff))/len(shuffled_beta_diff))
real_beta_diff = pd.DataFrame(real_beta_diff)
real_beta_diff["pval_beta_diff"] = real_beta_diff.b.apply(lambda x: get_pval_betadiff(x))
mild = real_beta_diff["pval_beta_diff"]<0.2
mild.index = d1.index
signif = real_beta_diff["pval_beta_diff"]<0.05
signif.index = d1.index
bonf = real_beta_diff["pval_beta_diff"]<0.05/d1.shape[0]
bonf.index = d1.index #not enough test..
#figure:
signif1 = (abs(d1.b - d2.b)>1) & (abs(d1.b)>abs(d2.b))
signif1.index = d1.index
signif2 = (abs(d1.b - d2.b)>1) & (abs(d1.b)<abs(d2.b))
signif2.index = d2.index
m = min(d1.b.min(),d2.b.min())
M = max(d1.b.max(),d2.b.max())
plt.figure(figsize=(4.5,4.5))
plt.scatter(x=d1.b, y=d2.b, color="tab:grey", edgecolor="black", alpha=0.75, label="Others")
plt.scatter(x=d1[signif1].b, y=d2[signif1].b, color="tab:orange", edgecolor="black", alpha=0.75, label="G->mRNA->Protein")
plt.scatter(x=d1[signif2].b, y=d2[signif2].b, color="tab:blue", edgecolor="black", alpha=0.75, label="G->Protein->mRNA")
plt.plot([m,M],[m,M], color="black", zorder=1, linestyle="--", linewidth=0.5)
plt.plot([m,M],[M,m], color="black", zorder=1, linestyle="--", linewidth=0.5)
plt.axvline(x=0, color="black", zorder=1, linestyle="--", linewidth=0.5)
plt.axhline(y=0, color="black", zorder=1, linestyle="--", linewidth=0.5)
plt.xlabel("Effect size\nexposure=mRNA, outcome=protein")
plt.ylabel("Effect size\nexposure=protein, outcome=mRNA")
plt.legend(title="Causal model:")
plt.xlim([m*1.05,M*1.05])
plt.ylim([m*1.05,M*1.05])
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/MR_effsizes.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/MR_effsizes.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#and the numbers:
ep = sum(signif1)
pe = sum(signif2)
oth = d1.shape[0]-pe-ep
plt.figure(figsize=(2,5))
plt.bar([0,1],[ep, pe], color=["tab:orange", "tab:blue"])
plt.text(x=0-.4, y=ep, s="n={0}".format(ep))
plt.text(x=1-.4, y=pe, s="n={0}".format(pe))
plt.xticks([0,1],["G->mRNA->Protein", "G->Protein->mRNA"], rotation=60)
plt.xlabel("Causal model")
plt.ylabel("Number of genes")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/MR_effsizes_N.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/MR_effsizes_N.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#and investigate the model1 and model2 a bit:
epgene = d1[signif1].gene_id
pegene = d1[signif2].gene_id
othgene = d1[~(signif1|signif2)].gene_id
#overlay their expression bin:
wb = pd.read_csv("~/Downloads/wbfrac_classification.tsv.gz", sep="\t", compression="gzip", index_col=0)
epgene = epgene.str.split("\\.").str[0]
pegene = pegene.str.split("\\.").str[0]
othgene = othgene.str.split("\\.").str[0]
wb["mr"] = "NA"
wb.loc[wb.index.intersection(epgene),"mr"] = "ep"
wb.loc[wb.index.intersection(pegene),"mr"] = "pe"
wb.loc[wb.index.intersection(othgene),"mr"] = "oth"
wb = wb[wb.mr!="NA"]
tb = wb.groupby(["mr", "wb_frac_bin"]).size().unstack().fillna(0).astype(int)
#Plot
tb = tb.loc[["ep","pe","oth"],["low", "mid", "high"]]
frac = (tb.T/tb.sum(axis=1)).T
err = np.sqrt((frac*(1-frac)).T / tb.sum(axis=1)).T
frac.columns.name = ""

frac.index = ["G->mRNA->Protein", "G->Protein->mRNA", "Others"]
palette = sns.cubehelix_palette(5, start=2)
frac.plot.bar(stacked=True, color=[palette[1], palette[2], palette[3]], rot=60, figsize=(4,5.5))
#adding errorbars:
plt.errorbar(np.arange(frac.shape[0]), frac.iloc[:,:-2].sum(axis=1), err.iloc[:,-2], fmt="none", color="black")
plt.errorbar(np.arange(frac.shape[0]), frac.iloc[:,:-1].sum(axis=1), err.iloc[:,-1], fmt="none", color="black")
plt.ylim([0,1])
plt.xlabel("Causal model")
plt.ylabel("Whole Blood mRNA\nexpression bin fraction")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/mr_wb_frac.png', dpi=500, bbox_inches='tight')
plt.savefig('/Users/qingbowang/Downloads/mr_wb_frac.pdf', dpi=500, bbox_inches='tight')
plt.clf()

#And finally, overlay the per-gene coloc classification
#cls = pd.read_csv('/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/es_ps_stats_pergene_updated.tsv',sep='\t')
cls = pd.read_csv('~/Downloads/es_ps_stats_pergene_updated.tsv',sep='\t')#in local
cls.index = cls.gene_id.str.split("\\.").str[0]
to_add = genes.str.split("\\.").str[0].difference(cls.index)
#Slow way to add rows:
dfadd = pd.DataFrame(0, index=to_add, columns = cls.columns)
cls = pd.concat([cls, dfadd], axis=0)
cls["mr"] = "NA"
cls.loc[epgene,"mr"] = "ep"
cls.loc[pegene,"mr"] = "pe"
cls.loc[othgene,"mr"] = "oth"
cls = cls[cls.mr!="NA"]
cls["coloc_evidence"] = 0
cls.loc[cls.coloc_strict+cls.coloc_lenient+cls.coloc_super_lenient>0,"coloc_evidence"] = 1
cls.loc[cls.coloc_strict+cls.coloc_lenient>0,"coloc_evidence"] = 2
tb = cls.groupby(["mr", "coloc_evidence"]).size().unstack().fillna(0).astype(int)
cls["es_evidence"] = cls.mrna>0
tb2 = cls.groupby(["mr", "es_evidence"]).size().unstack().fillna(0).astype(int)
cls["ps_evidence"] = cls.protein>0
tb3 = cls.groupby(["mr", "ps_evidence"]).size().unstack().fillna(0).astype(int)
cls["coloc_evidence"] = 0
cls.loc[cls.coloc_strict+cls.coloc_lenient+cls.coloc_super_lenient>0,"coloc_evidence"] = 1
cls.loc[cls.coloc_strict+cls.coloc_lenient>0,"coloc_evidence"] = 2
tb = cls.groupby(["mr", "coloc_evidence"]).size().unstack().fillna(0).astype(int)
tb = tb.loc[["ep","pe","oth"],[0,1,2]]
tb.columns = ["<0.01", "[0.01,0.1)", "0.1$\leq$"]
tb.index = ["G->mRNA->Protein", "G->Protein->mRNA", "Others"]
frac = (tb.T/tb.sum(axis=1)).T
err = np.sqrt((frac*(1-frac)).T / tb.sum(axis=1)).T
frac.columns.name = ""
frac.plot.bar(stacked=True, color=["tab:grey", "tab:olive", "tab:orange"], rot=60, figsize=(4,5.5))
#adding errorbars:
plt.errorbar(np.arange(frac.shape[0]), frac.iloc[:,:-2].sum(axis=1), err.iloc[:,-2], fmt="none", color="black")
plt.errorbar(np.arange(frac.shape[0]), frac.iloc[:,:-1].sum(axis=1), err.iloc[:,-1], fmt="none", color="black")
plt.ylim([0,1])
plt.xlabel("Causal model")
plt.ylabel("Colocalization evidence in CLPP")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/mr_clpp_frac.png', dpi=500, bbox_inches='tight')
plt.savefig('/Users/qingbowang/Downloads/mr_clpp_frac.pdf', dpi=500, bbox_inches='tight')
plt.clf()