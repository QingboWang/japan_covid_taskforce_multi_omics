#simply replicate in eQTLgen:
import pandas as pd
import time as tm
#our cis-eQTL fine-mapping results, based on hg19:
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=" ")
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep="\t")
e_fm["variant_id_hg19"] = e_fm.rsid.str.split("_").str[1]
e_fm["gene_id"] = e_fm.rsid.str.split("_").str[-1]
e_fm.set_index(["variant_id_hg19","gene_id"], inplace=True)
e_sus["variant_id_hg19"] = e_sus.rsid.str.split("_").str[1]
e_sus["gene_id"] = e_sus.rsid.str.split("_").str[-1]
e_sus.set_index(["variant_id_hg19","gene_id"], inplace=True)
ej = e_fm.join(e_sus[["rsid","pip"]], how="outer", rsuffix="_sus")
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej["vg"] = ej.rsid #rsidがsusie, fmどっち由来でもkeepする必要ありなので..
ej.loc[ej.rsid==0,"vg"] = ej.loc[ej.rsid==0,"rsid_sus"]
del ej["rsid"]
del ej["rsid_sus"]
ej.columns = ["pip_fm", "pip_sus", "vg"]
ej.reset_index(inplace=True)
ej["gene_id_unv"] = ej.gene_id.str.split("\\.").str[0]
ej.set_index(["variant_id_hg19", "gene_id_unv"], inplace=True)
#eQTLgen:
#full = 127341799 lines....
#First filter to fdr005 in their dataset:
gen = pd.read_csv("~/Desktop/resources/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz", sep='\t')
#gen = pd.read_csv("~/Desktop/resources/eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz", sep='\t')
print ("done reading eQTLgen info {0}".format(tm.ctime()))
gen["variant_hg19"] = gen.SNPChr.astype(str)+":"+gen.SNPPos.astype(str)+":"+gen.OtherAllele+":"+gen.AssessedAllele
gen["variant_hg19_flipped"] = gen.SNPChr.astype(str)+":"+gen.SNPPos.astype(str)+":"+gen.AssessedAllele+":"+gen.OtherAllele
gen.set_index(["variant_hg19", "Gene"], inplace=True)
tojoin = gen.Zscore.copy(deep=True)
tojoin.index.names = ej.index.names
ej = ej.join(tojoin, how="left")
gen.reset_index(inplace=True)
gen.set_index(["variant_hg19_flipped", "Gene"], inplace=True)
tojoin = gen.Zscore.copy(deep=True)
tojoin.index.names = ej.index.names
ej = ej.join(tojoin, how="left", rsuffix="_needsflip")
print ("done adding eQTLgen info {0}".format(tm.ctime()))
#add our JCTF z score:
ej.reset_index(inplace=True)
ej.set_index("vg", inplace=True)
df0 = []
for chr in list(range(1,23))+["X"]:
    df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t')
    df.index = df.variant_id + "_" + df.gene_id
    df = df.loc[df.index.intersection(ej.index),:]
    df["zscore_jctf"] = df.slope/df.slope_se
    df0.append(df.zscore_jctf)
    print ("done {0}, {1}".format(chr, tm.ctime()))
df0 = pd.concat(df0)
df0.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/eqtl_finemapORsusie_pip0001_zscore.tsv", sep="\t")
ej = ej.join(df0, how="left")
ej.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/eqtl_finemapORsusie_pip0001_foreQTLgencomparison.tsv", sep="\t")
#Started running 17:43, done 18:05 ish -> re-doing
#Done, working on this sumstats:
ej = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/eqtl_finemapORsusie_pip0001_foreQTLgencomparison.tsv", sep="\t", index_col=0)
#Zscore flip taiou
ej["zscore_eqtlgen"] = ej.Zscore
ej.loc[~ej.Zscore_needsflip.isna(),"zscore_eqtlgen"] = ej.loc[~ej.Zscore_needsflip.isna(),"Zscore_needsflip"]*-1
#e.g. simple replication rate:
plt.scatter(ej.zscore_jctf, ej.zscore_eqtlgen, alpha=0.5)
plt.show()
#Plot, per PIP
def bin_pip(pip):
    if pip<0.001: return 0
    elif pip<0.01: return 0.001
    elif pip<0.1: return 0.01
    elif pip<0.5: return 0.1
    elif pip<0.9: return 0.5
    elif pip<1: return 0.9
    elif pip==1: return 1
    else: return np.nan
ej["min_pip"] = np.minimum(ej.pip_fm, ej.pip_sus)
ej["min_pip_bin"] = ej.min_pip.apply(lambda x: bin_pip(x))
ej.min_pip_bin.value_counts()
colors = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
i = 0
for bin in ej.min_pip_bin.value_counts().sort_index().index:
    plt.scatter(ej[ej.min_pip_bin==bin].zscore_jctf, ej[ej.min_pip_bin==bin].zscore_eqtlgen, color=colors[i])
    i += 1
plt.show()

#First, plot overall density:
#taken from the last post of https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn

def density_scatter( x , y, ax = None, sort = True, bins = 20, **kwargs )   :
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins = bins, density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)
    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0
    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    ax.scatter( x, y, c=z, **kwargs )
    return ax
titles = ["[0.0001, 0.001)", "[0.001, 0.01)", "[0.01, 0.1)", "[0.1,0.5)", "[0.5,0.9)", "[0.9,1)", "1"]
ej_forplot = ej[~ej.zscore_eqtlgen.isna()]#removing NAs from the beginning
fig, axes = plt.subplots(2, 4, figsize=(8, 4.8), sharex=True, sharey=True)
x = ej_forplot.zscore_jctf
y = ej_forplot.zscore_eqtlgen
density_scatter(x, y, ax=axes[0, 0])
axes[0][0].axvline(x=0, linestyle="--", linewidth=1, color="black", zorder=-2)
axes[0][0].axhline(y=0, linestyle="--", linewidth=1, color="black", zorder=-2)
axes[0][0].text(x=-100, y=150, s="r={0}".format(np.round(stats.pearsonr(x, y)[0], 4)))
axes[0][0].set_title("Overall density")
cnt = 0
for bin in ej_forplot.min_pip_bin.value_counts().sort_index().index:
    if cnt<3:
        i=0
        j=cnt+1
    else:
        i=1
        j=cnt+1-4
    x = ej_forplot[ej_forplot.min_pip_bin==bin].zscore_jctf
    y = ej_forplot[ej_forplot.min_pip_bin==bin].zscore_eqtlgen
    axes[i][j].scatter(x, y, color=colors[cnt], alpha=0.7, edgecolor=colors[cnt])
    axes[i][j].axvline(x=0, linestyle="--", linewidth=1, color="black", zorder=-2)
    axes[i][j].axhline(y=0, linestyle="--", linewidth=1, color="black", zorder=-2)
    axes[i][j].text(x=-100, y=150, s="r={0}".format(np.round(stats.pearsonr(x,y)[0], 4)))
    axes[i][j].set_title("PIP: {0}".format(titles[cnt]), fontsize=12)
    cnt += 1
axes[1][1].set_xlabel("                                        Z score in JCTF")
axes[1][0].set_ylabel("                                                Z score in eQTLgen")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/vs_eQTLgen_perPIP.png", dpi=500)
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/vs_eQTLgen_perPIP.pdf", dpi=500)
plt.clf()
# also per PIP bin, lead var enrichment:
gen = pd.read_csv("~/Desktop/resources/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz", sep='\t')
gen["variant_hg19"] = gen.SNPChr.astype(str)+":"+gen.SNPPos.astype(str)+":"+gen.OtherAllele+":"+gen.AssessedAllele
gen["variant_hg19_flipped"] = gen.SNPChr.astype(str)+":"+gen.SNPPos.astype(str)+":"+gen.AssessedAllele+":"+gen.OtherAllele
lead = gen.sort_values(by=["Gene", "Pvalue"], ascending=True).groupby("Gene").head(1)
lead.set_index(["variant_hg19", "Gene"], inplace=True)
lead["is_lead"] = True
ej.set_index(["variant_id_hg19", "gene_id_unv"], inplace=True)
lead.index.names = ej.index.names
ej = ej.join(lead.is_lead, how="left")#no need to assume flip in this case. Too rare event.
ej["is_lead"] = ej["is_lead"].fillna(False)
ej["is_lead"].value_counts()

st = ej.groupby("min_pip_bin").is_lead.agg(["mean", "sem"])
st.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/eqtlgen_repr_fraclead.tsv", sep="\t")

colors = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
titles = ["[0.0001, 0.001)", "[0.001, 0.01)", "[0.01, 0.1)", "[0.1,0.5)", "[0.5,0.9)", "[0.9,1)", "1"]
plt.figure(figsize=(5,3.5))
for i in np.arange(st.shape[0]):
    plt.errorbar(i, st["mean"].values[i]*100, st["sem"].values[i]*100, fmt="o", c=colors[i])
plt.xticks(np.arange(st.shape[0]), titles, rotation=45)
plt.xlabel("PIP bin in JCTF")
plt.ylabel("%lead-variant in eQTLgen")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/vs_eQTLgen_leadvarenr.png", dpi=500)
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/vs_eQTLgen_leadvarenr.pdf", dpi=500)
plt.clf()



#Then the eQTLgen replication of our pQTLs that do not replicate in our eQTL:

ej = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/eqtl_finemapORsusie_pip0001_foreQTLgencomparison.tsv", sep="\t", index_col=0)
#Zscore flip cases
ej["zscore_eqtlgen"] = ej.Zscore
ej.loc[~ej.Zscore_needsflip.isna(),"zscore_eqtlgen"] = ej.loc[~ej.Zscore_needsflip.isna(),"Zscore_needsflip"]*-1
#And pQTL results as well:
allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eps_n998_pip0001.tsv.gz", sep="\t", index_col=[0,1])
#Focusing on the ones with pQTL PIP>0.001.
ej["variant_id_hg38"] = ej.index.str.split("_").str[0]
ej.set_index(["variant_id_hg38", "gene_id"], inplace=True)
allj = allj[allj.p_min_pip>0.001]
allj = allj.join(ej.zscore_eqtlgen, how="left")
allj["clpp"] = allj.e_min_pip*allj.p_min_pip
allj.loc[allj.clpp<0.001,"clpp"] = 0 #Clipping

#Test z score distribution:
allj["z_bin"] = 0
allj.loc[abs(allj.zscore_eqtlgen)>5,"z_bin"] = 5
allj.loc[abs(allj.zscore_eqtlgen)>10,"z_bin"] = 10
allj.loc[abs(allj.zscore_eqtlgen)>20,"z_bin"] = 20
allj.loc[abs(allj.zscore_eqtlgen)>50,"z_bin"] = 50
allj["p_bin"] = 0
allj.loc[allj.p_min_pip>0.001,"p_bin"] = 0.001
allj.loc[allj.p_min_pip>0.01,"p_bin"] = 0.01
allj.loc[allj.p_min_pip>0.1,"p_bin"] = 0.1
allj.loc[allj.p_min_pip>0.9,"p_bin"] = 0.9
allj["e_bin"] = 0
allj.loc[allj.e_min_pip>0.001,"e_bin"] = 0.001
allj.loc[allj.e_min_pip>0.01,"e_bin"] = 0.01
allj.loc[allj.e_min_pip>0.1,"e_bin"] = 0.1
allj.loc[allj.e_min_pip>0.9,"e_bin"] = 0.9
st = allj.groupby(["p_bin", "e_bin", "z_bin"]).size().unstack().fillna(0).astype(int)

#CLPP style is easier
allj["clpp_bin"] = 0
allj.loc[allj.clpp>0.001,"clpp_bin"] = 0.001
allj.loc[allj.clpp>0.01,"clpp_bin"] = 0.01
allj.loc[allj.clpp>0.1,"clpp_bin"] = 0.1
allj.loc[allj.clpp>0.9,"clpp_bin"] = 0.9
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib.colors import LogNorm
import seaborn as sns
st = allj.groupby(["e_bin", "z_bin"]).size().unstack().fillna(0).astype(int)
#Plot the raw number:
st.index = ["<0.001", "[0.001,0.01)", "[0.01,0.1)", "[0.1,0.9)", "0.9$\leq$"]
st.columns = ["<5\n(or missing)", "[5,10)","[10,20)", "[20,50)", "50$\leq$"]
tb = st.T.iloc[::-1,:] #to match the order
log_norm = LogNorm()
plt.figure(figsize=(6,6))
sns.heatmap(tb+1, annot=tb, fmt="d", square=True, linewidths=.5, norm=log_norm,
            cmap="viridis", cbar_kws={'label': 'count',
                                      "shrink": .6})
plt.yticks(rotation=35, fontsize=13)
plt.xticks(rotation=35, fontsize=13)
plt.xlabel("eQTL PIP bin, our data", fontsize=15)
plt.ylabel("|z| in eQTLgen", fontsize=15)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/eqtlgen_pip_heatmap.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/eqtlgen_pip_heatmap.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#And column-normalized percentage:
stnorm = np.round((st.T/st.sum(axis=1)).T * 100, 1)
tb = stnorm.T.iloc[::-1,:] #to match the order
plt.figure(figsize=(6,6))
sns.heatmap(tb, annot=tb, fmt=".1f", square=True, linewidths=.5,
            cmap="viridis", cbar_kws={'label': 'Column-normalized %',
                                      "shrink": .6})
plt.yticks(rotation=35, fontsize=13)
plt.xticks(rotation=35, fontsize=13)
plt.xlabel("eQTL PIP bin, our data", fontsize=15)
plt.ylabel("|z| in eQTLgen", fontsize=15)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/eqtlgen_pip_heatmap_norm.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/eqtlgen_pip_heatmap_norm.pdf', bbox_inches='tight', dpi=500)
plt.clf()

allj["pip_eqtlgen_approx"] = 0
allj.loc[allj.zscore_eqtlgen>5,"pip_eqtlgen_approx"] = 0.5
allj.loc[allj.zscore_eqtlgen>10,"pip_eqtlgen_approx"] = 1 #Let's be super lenient as a test case
allj["clpp_whatif"] = allj.p_min_pip * np.maximum(allj.e_min_pip, allj.pip_eqtlgen_approx)
allj["clpp_whatif_bin"] = 0
allj.loc[allj.clpp_whatif>0.001,"clpp_whatif_bin"] = 0.001
allj.loc[allj.clpp_whatif>0.01,"clpp_whatif_bin"] = 0.01
allj.loc[allj.clpp_whatif>0.1,"clpp_whatif_bin"] = 0.1
allj.loc[allj.clpp_whatif>0.9,"clpp_whatif_bin"] = 0.9
st = allj.groupby(["clpp_bin", "clpp_whatif_bin"]).size().unstack().fillna(0).astype(int)
#Plot the raw number:
st.index = ["<0.001", "[0.001,0.01)", "[0.01,0.1)", "[0.1,0.9)", "0.9$\leq$"]
st.columns = st.index
tb = st.T.iloc[::-1,:] #to match the order
log_norm = LogNorm()
plt.figure(figsize=(6,6))
sns.heatmap(tb+1, annot=tb, fmt="d", square=True, linewidths=.5, norm=log_norm,
            cmap="viridis", cbar_kws={'label': 'count',
                                      "shrink": .6})
plt.yticks(rotation=35, fontsize=13)
plt.xticks(rotation=35, fontsize=13)
plt.xlabel("CLPP bin, our data", fontsize=15)
plt.ylabel("CLPP bin, if eQTLgen-\nsignificant variants are PIP~=1", fontsize=15)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/eqtlgen_whatif_heatmap.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/eqtlgen_whatif_heatmap.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#And column-normalized percentage:
stnorm = np.round((st.T/st.sum(axis=1)).T * 100, 1)
tb = stnorm.T.iloc[::-1,:] #to match the order
plt.figure(figsize=(6,6))
sns.heatmap(tb, annot=tb, fmt=".1f", square=True, linewidths=.5,
            cmap="viridis", cbar_kws={'label': 'Column-normalized %',
                                      "shrink": .6})
plt.yticks(rotation=35, fontsize=13)
plt.xticks(rotation=35, fontsize=13)
plt.xlabel("CLPP bin, our data", fontsize=15)
plt.ylabel("CLPP bin, if eQTLgen-\nsignificant variants are PIP~=1", fontsize=15)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/eqtlgen_whatif_heatmap_norm.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/eqtlgen_whatif_heatmap_norm.pdf', bbox_inches='tight', dpi=500)
plt.clf()



