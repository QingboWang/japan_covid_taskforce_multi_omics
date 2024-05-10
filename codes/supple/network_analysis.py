#Analyses based on https://github.com/mortazavilab/PyWGCNA#readme

#load the df
import pandas as pd
import PyWGCNA
from matplotlib import pyplot as plt
#write mrna and protein separately, and also full versus intersection:
fn="~/Desktop/taskforce_n1102/n1005.protein.invnormed.expression.bed"
pe = pd.read_csv(fn, sep='\t', index_col=[0,1,2,3])
fn="~/Desktop/taskforce_n1102/n1019.expression.bed.gz"
me = pd.read_csv(fn, sep='\t', index_col=[0,1,2,3])

me.index = me.index.get_level_values(3)
pe.index = pe.index.get_level_values(3)
me.columns.names = ["sample_id"]
pe.columns.names = ["sample_id"]
#full:
me.T.to_csv("~/Desktop/taskforce_n1102/network/mrna_n1005_fornetwork.csv", sep=",")
pe.T.to_csv("~/Desktop/taskforce_n1102/network/protein_n1005_fornetwork.csv", sep=",")
#filtered to intersection:
me = me.loc[pe.index.intersection(me.index), pe.columns]
pe = pe.loc[me.index, :]
me.index = me.index.str.split("\\.").str[0]
pe.index = pe.index.str.split("\\.").str[0]
me.T.to_csv("~/Desktop/taskforce_n1102/network/mrna_n1005_intersection_fornetwork.csv", sep=",")
pe.T.to_csv("~/Desktop/taskforce_n1102/network/protein_n1005_intersection_fornetwork.csv", sep=",")

#WGCNA; focusing on 2000 genes should be good.
geneExp = '/Users/qingbowang/Desktop/taskforce_n1102/network/mrna_n1005_intersection_fornetwork.csv'
obj_e = PyWGCNA.WGCNA(name='jctf_rna', species='human',
                              geneExpPath=geneExp, TPMcutoff = -100,
                              save=True, outputPath="/Users/qingbowang/Desktop/taskforce_n1102/network/")
obj_e.geneExpr.to_df().head(5)
obj_e.preprocess() #did this go well? Not sure but for now, proceed...
#find modules
obj_e.findModules()
#just save and check things:
obj_e.saveWGCNA()

#protein:
geneExp = '/Users/qingbowang/Desktop/taskforce_n1102/network/protein_n1005_intersection_fornetwork.csv'
obj_p = PyWGCNA.WGCNA(name='jctf_protein', species='human',
                              geneExpPath=geneExp, TPMcutoff = -100,
                              save=True, outputPath="/Users/qingbowang/Desktop/taskforce_n1102/network/")
obj_p.geneExpr.to_df().head(5)
obj_p.preprocess() #did this go well? Not sure but for now, proceed...
#find modules
obj_p.findModules()
#just save and check things:
obj_p.saveWGCNA()

#comparison:
obj_e = PyWGCNA.readWGCNA("/Users/qingbowang/Desktop/taskforce_n1102/network/jctf_rna.p")
obj_p = PyWGCNA.readWGCNA("/Users/qingbowang/Desktop/taskforce_n1102/network/jctf_protein.p")

comparison = PyWGCNA.compareNetworks(PyWGCNAs = [obj_e, obj_p])
#e.g. similarity:
comparison.jaccard_similarity.head(5)
#visualize:
color = {"jctf_rna": "lightblue",
         "jctf_protein": "lightgreen"}
comparison.plotJaccardSimilarity(color=color,
                                 cutoff=0.05,
                                 plot_format="pdf",
                                 file_name="/Users/qingbowang/Desktop/taskforce_n1102/network/figures/jaccard_similarity")
comparison.plotHeatmapComparison(plot_format="pdf",
                                 file_name="/Users/qingbowang/Desktop/taskforce_n1102/network/figures/heatmap_similarity")
comparison.plotBubbleComparison(color=color,
                                figsize=(10,10),
                                plot_format="pdf",
                                file_name="/Users/qingbowang/Desktop/taskforce_n1102/network/figures/bubble_similarity")
#Also png:
comparison.plotJaccardSimilarity(color=color,
                                 cutoff=0.05,
                                 plot_format="png",
                                 file_name="/Users/qingbowang/Desktop/taskforce_n1102/network/figures/jaccard_similarity")
comparison.plotHeatmapComparison(plot_format="png",
                                 file_name="/Users/qingbowang/Desktop/taskforce_n1102/network/figures/heatmap_similarity")
comparison.plotBubbleComparison(color=color,
                                figsize=(10,10),
                                plot_format="png",
                                file_name="/Users/qingbowang/Desktop/taskforce_n1102/network/figures/bubble_similarity")

#GO term:
#(“snow”, "darkgrey" and “dimgrey” in mRNA, “black” and “darkgrey” in protein expression)
geneList = PyWGCNA.getGeneList(dataset='hsapiens_gene_ensembl',
                               attributes=['ensembl_gene_id',
                                           'external_gene_name',
                                           'gene_biotype'],
                               maps=['gene_id', 'gene_name', 'gene_biotype']) #taking some time....
gene_set_library = ["Reactome_2022"]
obj_e.figureType = "png"
obj_p.figureType = "png"
obj_e.updateGeneInfo(geneList)
obj_p.updateGeneInfo(geneList)

obj_e.functional_enrichment_analysis(type="GO",
                                             moduleName="darkgrey",
                                             sets=gene_set_library,
                                             p_value=0.05,
                                             file_name="reac_mrna_darkgrey")
obj_e.functional_enrichment_analysis(type="GO",
                                             moduleName="dimgrey",
                                             sets=gene_set_library,
                                             p_value=0.05,
                                             file_name="reac_mrna_dimgrey")
obj_e.functional_enrichment_analysis(type="GO",
                                             moduleName="snow",
                                             sets=gene_set_library,
                                             p_value=0.05,
                                             file_name="reac_mrna_snow")
obj_p.functional_enrichment_analysis(type="GO",
                                             moduleName="darkgrey",
                                             sets=gene_set_library,
                                             p_value=0.05,
                                             file_name="reac_protein_darkgrey")
obj_p.functional_enrichment_analysis(type="GO",
                                             moduleName="black",
                                             sets=gene_set_library,
                                             p_value=0.05,
                                             file_name="reac_protein_black")

#Neg control:
obj_e.functional_enrichment_analysis(type="GO",
                                             moduleName="white",
                                             sets=gene_set_library,
                                             p_value=0.05,
                                             file_name="reac_mrna_white")
#type="REACTOME" somehow gets null... fine..

#Scatter between modules:
#For that, p=1まで撮り直す;
obj_e.functional_enrichment_analysis(type="GO",
                                             moduleName="darkgrey",
                                             sets=gene_set_library,
                                             p_value=1,
                                             file_name="reac_mrna_darkgrey_full")
obj_e.functional_enrichment_analysis(type="GO",
                                             moduleName="dimgrey",
                                             sets=gene_set_library,
                                             p_value=1,
                                             file_name="reac_mrna_dimgrey_full")
obj_e.functional_enrichment_analysis(type="GO",
                                             moduleName="snow",
                                             sets=gene_set_library,
                                             p_value=1,
                                             file_name="reac_mrna_snow_full")
obj_p.functional_enrichment_analysis(type="GO",
                                             moduleName="darkgrey",
                                             sets=gene_set_library,
                                             p_value=1,
                                             file_name="reac_protein_darkgrey_full")
obj_p.functional_enrichment_analysis(type="GO",
                                             moduleName="black",
                                             sets=gene_set_library,
                                             p_value=1,
                                             file_name="reac_protein_black_full")
#Neg control:
obj_e.functional_enrichment_analysis(type="GO",
                                             moduleName="white",
                                             sets=gene_set_library,
                                             p_value=1,
                                             file_name="reac_mrna_white_full")
e_darkgrey = pd.read_csv("~/Desktop/taskforce_n1102/network/figures/GO/reac_mrna_darkgrey_full/Reactome_2022.human.enrichr.reports.txt", sep="\t", index_col=1)
e_dimgrey = pd.read_csv("~/Desktop/taskforce_n1102/network/figures/GO/reac_mrna_dimgrey_full/Reactome_2022.human.enrichr.reports.txt", sep="\t", index_col=1)
e_snow = pd.read_csv("~/Desktop/taskforce_n1102/network/figures/GO/reac_mrna_snow_full/Reactome_2022.human.enrichr.reports.txt", sep="\t", index_col=1)
e_white = pd.read_csv("~/Desktop/taskforce_n1102/network/figures/GO/reac_mrna_white_full/Reactome_2022.human.enrichr.reports.txt", sep="\t", index_col=1)
p_darkgrey = pd.read_csv("~/Desktop/taskforce_n1102/network/figures/GO/reac_protein_darkgrey_full/Reactome_2022.human.enrichr.reports.txt", sep="\t", index_col=1)
p_black = pd.read_csv("~/Desktop/taskforce_n1102/network/figures/GO/reac_protein_black_full/Reactome_2022.human.enrichr.reports.txt", sep="\t", index_col=1)
len(e_darkgrey.index.intersection(e_dimgrey.index)) #Why still limited intersection? Anyways
uni = pd.concat([e_darkgrey["P-value"], e_dimgrey["P-value"], e_snow["P-value"], e_white["P-value"], p_darkgrey["P-value"], p_black["P-value"]], axis=1)
uni = uni.fillna(1)
uni.columns = ["mRNA:darkgrey", "mRNA:dimgrey", "mRNA:snow", "mRNA:white (negative)", "protein:darkgrey", "protein:black"]
uni = np.log10(uni)*-1
uni.to_csv("~/Desktop/taskforce_n1102/network/figures/GO/pval_comparison.tsv", sep="\t")


#e.g.
plt.scatter(uni.iloc[:,0], uni.iloc[:,-1])
plt.show()
#Plot overall, in local
from matplotlib import pyplot as plt
import pandas as pd
import time as tm
import numpy as np
from scipy import stats
import matplotlib.colors as mcolors
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
uni = pd.read_csv("~/Downloads/pval_comparison.tsv", sep="\t", index_col=0)
uni.columns = uni.columns.str.replace("\n", " ") #改行いらんかったわ.
for i, col1 in enumerate(uni.columns):
    for j, col2 in enumerate(uni.columns):
        plt.figure(figsize=(4,4))
        plt.xlabel(col2 + " REACTOME -log10(p)")
        plt.ylabel(col1 + " REACTOME -log10(p)")
        x = uni[col2]
        y = uni[col1]
        if ("white" in col1) | ("white" in col2):
            plt.scatter(x, y, alpha=0.75, edgecolor="tab:grey", color="tab:grey")
        else:
            plt.scatter(x, y, alpha=0.75, edgecolor="tab:blue", color="tab:blue")
        plt.axvline(x=0, color="tab:grey", linestyle="--", linewidth=0.5, zorder=-2)
        plt.axhline(y=0, color="tab:grey", linestyle="--", linewidth=0.5, zorder=-2)
        plt.plot([0, x.max()*1.1], [0, y.max()*1.1], color="tab:grey", linestyle="--", linewidth=0.5, zorder=-2)
        nonna = (~np.isnan(x)) & (~np.isnan(y))
        r, p = stats.pearsonr(x[nonna], y[nonna])
        plt.text(x=0, y=y.max(), s="r={0}".format(np.round(r, 4)), ha="left", va="top")
        plt.tight_layout()
        plt.savefig("/Users/qingbowang/Downloads/go_plots/x_{0}_y_{1}.png".format(col2, col1), dpi=500)
        plt.savefig("/Users/qingbowang/Downloads/go_plots/x_{0}_y_{1}.pdf".format(col2, col1), dpi=500)
        plt.clf()
dfe = obj_e.datExpr.var
dfp = obj_p.datExpr.var
#overlay the WB expression bin:
wb = pd.read_csv("~/Desktop/taskforce_n1102/n1300/wbfrac_classification.tsv.gz", sep="\t", compression="gzip", index_col=0)
dfe = dfe.join(wb, how="left")
dfp = dfp.join(wb, how="left")
dfe.to_csv("~/Desktop/taskforce_n1102/network/wgcna_df_e.tsv", sep="\t")
dfp.to_csv("~/Desktop/taskforce_n1102/network/wgcna_df_p.tsv", sep="\t")
ste = dfe.groupby(["moduleColors", "wb_frac_bin"]).size().unstack().fillna(0).astype(int)
stp = dfp.groupby(["moduleColors", "wb_frac_bin"]).size().unstack().fillna(0).astype(int)
ste = ste.loc[:,["high", "mid", "low"]]
stp = stp.loc[:,["high", "mid", "low"]]
(ste.T/ste.sum(axis=1)).T

#The not-too-different distribution could be interesting:
plte = ste.loc[["snow", "darkgrey", "dimgrey"],:]
pltp = stp.loc[["black", "darkgrey"],:]
plte.index = "mRNA:" + plte.index
pltp.index = "protein:" + pltp.index
#background:
pltbg = dfe.wb_frac_bin.value_counts()
pltbg.name = "Background"

toplot = pd.concat([plte.T, pltp.T, pltbg], axis=1).T
toplot.to_csv("~/Desktop/taskforce_n1102/network/wgcna_wbexp_toplot.tsv", sep="\t")
N = toplot.sum(axis=1)
for idx in toplot.index:
    toplot.rename(index={idx:idx+"\nn={0}".format(N[idx])}, inplace=True)
N = toplot.sum(axis=1) #re-naming
frac = (toplot.T/N).T
err = (np.sqrt( (frac*(1-frac)).T / N)).T
palette = sns.cubehelix_palette(5, start=2)

frac.loc[:,::-1].plot.bar(stacked=True, color=[palette[1], palette[2], palette[3]], rot=60, figsize=(4.5,4.5))
#adding errorbars:
plt.errorbar(np.arange(frac.shape[0]), frac.iloc[:,-2:].sum(axis=1), err.iloc[:,-2], fmt="none", color="black")
plt.errorbar(np.arange(frac.shape[0]), frac.iloc[:,-1:].sum(axis=1), err.iloc[:,-1], fmt="none", color="black")
plt.ylim([0,1])
plt.xlabel("Module")
plt.ylabel("Whole Blood mRNA\nexpression bin fraction")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/network/wbfrac_permodule.png', dpi=500, bbox_inches='tight')
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/network/wbfrac_permodule.pdf', dpi=500, bbox_inches='tight')
plt.clf()