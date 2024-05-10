import pandas as pd
import numpy as np
import time as tm
import glob

def run_my_susie_coloc(gn="ENSG00000175164.15"):
    print ("starting {0}, {1}".format(gn, tm.ctime()))
    cs = pd.read_csv("~/Desktop/taskforce_n1102/n1300/fmoutputs/eqtl_n998/{0}/{0}_e_susie_cs.txt".format(gn), sep="\t")
    alpha = pd.read_csv("~/Desktop/taskforce_n1102/n1300/fmoutputs/eqtl_n998/{0}/{0}_e_susie_alpha.tsv".format(gn),sep="\t")
    cs.set_index("cs", inplace=True)
    cs.sort_index(inplace=True)
    bf_df_e = []  # dataframe of BFs
    for csnum in cs.index:
        cols = str(cs.loc[csnum, "variable"]).split(",")
        cols = ["V" + sub for sub in cols]
        pips = alpha.loc[csnum, cols]
        log10bf = cs.loc[csnum, "cs_log10bf"]
        if np.isinf(log10bf):
            log10bf = 308  # to replace inf with 10**308
            print ("Inf in the eQTL log10(bf) for csID {0} detected, forcing to 1e308".format(csnum))
        sumbf = 10 ** log10bf
        bfs = sumbf * pips #assuming sum(alpha) in 95%CS ~= 1
        bf_df_e.append(bfs)
    bf_df_e = pd.concat(bf_df_e, axis=1).fillna(0)
    cs = pd.read_csv("~/Desktop/taskforce_n1102/n1300/fmoutputs/pqtl_n998/{0}/{0}_p_susie_cs.txt".format(gn), sep="\t")
    alpha = pd.read_csv("~/Desktop/taskforce_n1102/n1300/fmoutputs/pqtl_n998/{0}/{0}_p_susie_alpha.tsv".format(gn),sep="\t")
    cs.set_index("cs", inplace=True)
    cs.sort_index(inplace=True)
    bf_df_p = []  # dataframe of BFs
    for csnum in cs.index:
        cols = str(cs.loc[csnum, "variable"]).split(",")
        cols = ["V" + sub for sub in cols]
        pips = alpha.loc[csnum, cols]
        log10bf = cs.loc[csnum, "cs_log10bf"]
        if np.isinf(log10bf):
            log10bf = 308  # to replace inf with 10**308
            print("Inf in the pQTL log10(bf) for csID {0} detected, forcing to 1e308".format(csnum))
        sumbf = 10 ** log10bf
        bfs = sumbf * pips
        bf_df_p.append(bfs)
    bf_df_p = pd.concat(bf_df_p, axis=1).fillna(0)
    bf_df_e.columns = "e" + bf_df_e.columns.astype(str)
    bf_df_p.columns = "p" + bf_df_p.columns.astype(str)
    bf_df = pd.concat([bf_df_e, bf_df_p], axis=1).fillna(0)
    # and split to e and p again:
    bf_df_e = bf_df.loc[:, bf_df.columns.str.startswith("e")]
    bf_df_p = bf_df.loc[:, bf_df.columns.str.startswith("p")]
    # now calc. each PP:
    # 1 = e, 2 = pQTL
    p12 = 10 ** -6
    p1 = 10 ** -4
    p2 = 10 ** -4
    outdf = {}
    for col_e in bf_df_e.columns:
        for col_p in bf_df_p.columns:
            warnings = ""
            bf1 = bf_df_e.loc[:, col_e]
            bf2 = bf_df_p.loc[:, col_p]
            #pp0 = 1
            #pp1 = p1 * sum(bf1)
            #pp2 = p2 * sum(bf2) #These by definition would not exceed the original BF, so would not be inf
            #pp3 = p1 * p2 * (sum(np.outer(bf1, bf2)).sum() - sum(bf1 * bf2))
            #pp4 = p12 * sum(bf1 * bf2)
            pp0 = 1e-100
            pp1 = p1 * sum(bf1*1e-100)
            pp2 = p2 * sum(bf2*1e-100)
            diag = sum((bf1*1e-50) * (bf2*1e-50))
            if np.isposinf(diag): #if large enough to immediately call colocalization
                diag = 10 ** 308  # to replace max
                pp3 = 0 #force PP4=1 and pp3=0
                warnings += "pp4Inf"
            else:
                pp3 = p1 * p2 * (sum(np.outer((bf1*1e-50), (bf2*1e-50))).sum() - sum((bf1*1e-50) * (bf2*1e-50)))
            pp4 = p12 * diag
            if np.isinf(pp3):
                pp3 = p1 * p2 * 10**308
                warnings += "pp3Inf"
            print (pp0, pp1, pp2, pp3, pp4)
            sumpp = pp0 + pp1 + pp2 + pp3 + pp4
            pp0 = pp0 / sumpp
            pp1 = pp1 / sumpp
            pp2 = pp2 / sumpp
            pp3 = pp3 / sumpp
            pp4 = pp4 / sumpp
            outdf[(col_e, col_p)] = [pp0, pp1, pp2, pp3, pp4, warnings]
    outdf = pd.DataFrame(outdf).T
    outdf.index.names = ["e_cs", "p_cs"]
    outdf.columns = ["pp0","pp1", "pp2", "pp3", "pp4", "warnings"]
    print ("done {0}, {1}".format(gn, tm.ctime()))
    return (outdf)

genes = pd.read_csv("~/Desktop/taskforce_n1102/n1300/n998_mrna.protein_regressedout.bed", sep='\t', usecols=[0,1,2,3]).gene_id
cnt = 0
for gn in genes:
    try:
        cs = pd.read_csv("~/Desktop/taskforce_n1102/n1300/fmoutputs/eqtl_n998/{0}/{0}_e_susie_cs.txt".format(gn), sep="\t")
        csp = pd.read_csv("~/Desktop/taskforce_n1102/n1300/fmoutputs/pqtl_n998/{0}/{0}_p_susie_cs.txt".format(gn), sep="\t")
    except:
        print ("No CS for {0}, skipping".format(gn))
    else:
        try:
            sc = run_my_susie_coloc(gn)
            print (sc)
            sc.to_csv("~/Desktop/taskforce_n1102/n1300/susie_coloc_results/susie_coloc_{0}.tsv".format(gn), sep="\t")
        except:
            print ("Error for {0}, skipping".format(gn))
    if cnt%100==0:
        print ("done {0}th, {1}".format(cnt, tm.ctime()))
    cnt += 1
#parse the result for n=2211 genes:
genes = pd.read_csv("~/Desktop/taskforce_n1102/n1300/n998_mrna.protein_regressedout.bed", sep='\t', usecols=[0,1,2,3]).gene_id
pp3 = []
pp4 = []
for gn in genes:
    try:
        df = pd.read_csv("~/Desktop/taskforce_n1102/n1300/susie_coloc_results/susie_coloc_{0}.tsv".format(gn), sep="\t")
        pp3.append(df.pp3.max())
        pp4.append(df.pp4.max())
    except:
        pp3.append(-1)  # for originally not defined
        pp4.append(-1) #for originally not defined
pp4 = pd.Series(pp4)
pp4.index = genes
pp4.sort_values()
pp3 = pd.Series(pp3)
pp3.index = genes
pp3.sort_values()

pp4.to_csv("~/Desktop/taskforce_n1102/n1300/susie_coloc_results/all_pp4.tsv", sep="\t")


#Also simple ABF style coloc:
W = 0.15 #for ABF prior
genes = pd.read_csv("~/Desktop/taskforce_n1102/n1300/n998_mrna.protein_regressedout.bed", sep='\t', usecols=[0,1,2,3]).gene_id
for gn in genes:
    print ("starting {0}".format(gn))
    try:
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/eqtl_n998/{0}_e_fminput.z".format(gn),sep=" ", index_col=0)
        zp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gn),sep=" ", index_col=0)
        z["z"] = z.beta / z.se
        z["V"] = z.se ** 2
        z["r"] = (W / (W + z.V))
        z["log10abf_e"] = (np.log(np.sqrt(1 - z.r)) + (z.z ** 2 / 2 * z.r)) / np.log(10)  # maxが1e300になるよう平行移動して、underflowしているものはtrivialということで0に回す.
        zp["z"] = zp.beta / zp.se
        zp["V"] = zp.se ** 2
        zp["r"] = (W / (W + zp.V))
        zp["log10abf_p"] = (np.log(np.sqrt(1 - zp.r)) + (zp.z ** 2 / 2 * zp.r)) / np.log(10)  # maxが1e300になるよう平行移動して、underflowしているものはtrivialということで0に回す.
        log10abfs = pd.concat([z.log10abf_e,zp.log10abf_p], axis=1)
        log10abfs.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/abf_n998/{0}_log10abf.tsv".format(gn),sep="\t")
    except:
        print ("error for {0}".format(gn))
#and then the previous style coloc:
m = []
for gn in genes:
    print ("starting {0}".format(gn))
    abfs = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/abf_n998/{0}_log10abf.tsv".format(gn), sep="\t", index_col=0)
    m = m + [abfs.max().max()]
m = pd.Series(m)
m.hist(bins = 100)
plt.show()
outdf = {}
cnt = 0
for gn in genes:
    print ("starting {0}, {1}th, {2}".format(gn, cnt, tm.ctime()))
    try:
        abfs = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/abf_n998/{0}_log10abf.tsv".format(gn),sep="\t")
        bf1 = 10**abfs.log10abf_e.apply(lambda x: min(x, 300))
        bf2 = 10**abfs.log10abf_p.apply(lambda x: min(x, 300))
        p12 = 10 ** -6
        p1 = 10 ** -4
        p2 = 10 ** -4
        pp0 = 1e-100
        pp1 = p1 * sum(bf1 * 1e-100)
        pp2 = p2 * sum(bf2 * 1e-100)
        diag = sum((bf1 * 1e-50) * (bf2 * 1e-50))
        if np.isposinf(diag):  # if large enough to immediately call colocalization
            diag = 10 ** 308  # to replace max
            pp3 = 0  # force PP4=1 and pp3=0
        else:
            pp3 = p1 * p2 * (sum(np.outer((bf1 * 1e-50), (bf2 * 1e-50))).sum() - sum((bf1 * 1e-50) * (bf2 * 1e-50)))
        pp4 = p12 * diag
        if np.isinf(pp3):
            pp3 = p1 * p2 * 10 ** 308
        sumpp = pp0 + pp1 + pp2 + pp3 + pp4
        pp0 = pp0 / sumpp
        pp1 = pp1 / sumpp
        pp2 = pp2 / sumpp
        pp3 = pp3 / sumpp
        pp4 = pp4 / sumpp
        outdf[gn] = [pp0, pp1, pp2, pp3, pp4]
    except:
        print ("{0} error".format(gn))
    cnt += 1
outdf = pd.DataFrame(outdf)
outdf.index = ["pp0", "pp1", "pp2", "pp3", "pp4"]
outdf.T.to_csv("~/Desktop/taskforce_n1102/n1300/susie_coloc_results/simplecoloc_all_pp4.tsv", sep="\t")


#Now compare SuSiE coloc vs coloc vs our CLPP
allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eps_n998_pip001.tsv.gz", sep="\t")
allj["clpp"] = allj.e_min_pip*allj.p_min_pip
clpp = allj.groupby("gene_id").clpp.max()
simple_coloc = pd.read_csv("~/Desktop/taskforce_n1102/n1300/susie_coloc_results/simplecoloc_all_pp4.tsv", sep="\t", index_col=0)
susie_coloc = pd.read_csv("~/Desktop/taskforce_n1102/n1300/susie_coloc_results/all_pp4.tsv", sep="\t", index_col=0)
st = pd.concat([clpp, simple_coloc.pp4, susie_coloc.replace(-1,0)], axis=1)
st.columns = ["clpp", "simple_coloc", "susie_coloc"]
st = st[~st.simple_coloc.isna()]#simple coloc PP4 is na -> x chromsome. Remove from the bunbo
st = st.fillna(0)
#x = threshold (although acknowledging CLPP != PP.H4), y = %genes
stmean = {}
stsem = {}
#thresholds = [0.001]+list(np.arange(0.01,1.01, 0.01))
thresholds = np.logspace(-3,0,31) #log linear
for thres in thresholds:
    stmean[thres] = (st > thres).mean(axis=0)
    stsem[thres] = (st > thres).sem(axis=0)
stmean = pd.DataFrame(stmean)
stsem = pd.DataFrame(stsem)
#Plot:
rorder = ["susie_coloc", "simple_coloc"]
leg = ["SuSiE-coloc", "(simple) coloc"]
colors = ["tab:green", "tab:blue"]
x = stsem.columns
plt.figure(figsize=(6,4.2))
for i in range(3):
    plt.errorbar(x, stmean.loc[rorder[i],:], stsem.loc[rorder[i],:], color=colors[i], fmt="o", label=leg[i])
plt.xlim([0.00091,1.1])
plt.xscale("log")
plt.xticks(thresholds[::10], thresholds[::10], fontsize=14)
plt.xlabel("Threshold (CLPP or PP.H4)", fontsize=14)
plt.ylabel("Fraction of genes with colocalization\nsignal above threshold", fontsize=14)
plt.legend(title="Method:", loc="upper right")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/clpp_coloc_comparison.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/clpp_coloc_comparison.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#3x3 scatter plot with pearson r
from scipy import stats
import matplotlib.colors as mcolors
zmin = 0
zmax = 1
fig, axes = plt.subplots(3, 3, figsize=(6, 6), sharex=True, sharey=True)
for i, col1 in enumerate(st.columns):
    for j, col2 in enumerate(st.columns):
        if i == 2:
            axes[i, j].set_xlabel(col2)
        if j == 0:
            axes[i, j].set_ylabel(col1)
        if i>j:
            x = st[col2]
            y = st[col1]
            # Create density scatter plots
            hb = axes[i, j].hexbin(x, y, gridsize=20, cmap='binary', norm=mcolors.LogNorm())
        elif i<j:
            x = st[col2]
            y = st[col1]
            axes[i,j].scatter(x, y, alpha=0.75, edgecolor="tab:blue", color="tab:blue")
            axes[i,j].axvline(x=0, color="tab:grey", linestyle="--", linewidth=0.5, zorder=-2)
            axes[i,j].axhline(y=0, color="tab:grey", linestyle="--", linewidth=0.5, zorder=-2)
            axes[i,j].plot([zmin, zmax], [zmin, zmax], color="tab:grey", linestyle="--", linewidth=0.5, zorder=-2)
            nonna = (~np.isnan(x)) & (~np.isnan(y))
            r,p = stats.pearsonr(x[nonna], y[nonna])
            axes[i,j].text(x=zmax/2, y=zmax/2.2, s="r={0}".format(np.round(r, 4)), ha="left", va="top")
axes[0,0].set_xlim([zmin, zmax])
axes[0,0].set_ylim([zmin, zmax])
#fig.suptitle("Colocalization test comparisons")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Downloads/vscoloc_corr.png", dpi=500)
plt.savefig("/Users/qingbowang/Downloads/vscoloc_corr.pdf", dpi=500)
plt.clf()

#Bin is definitely easier to digest...
stbin = st.copy(deep=True)
stbin[-1<st] = 0
stbin[0.1<st] = 0.1
stbin[0.9<st] = 0.9
tb1 = stbin[["clpp","simple_coloc"]].value_counts().unstack().fillna(0).astype(int)
tb2 = stbin[["simple_coloc","susie_coloc"]].value_counts().unstack().fillna(0).astype(int)
tb3 = stbin[["susie_coloc","clpp"]].value_counts().unstack().fillna(0).astype(int)
#e.g. x=CLPP, y = SuSiE coloc
fracs = tb3/tb3.sum()
fracs.T.plot(kind='bar', stacked=True, figsize=(8, 6))
plt.show()
#そうする;
zmin = 0
zmax = 1
from matplotlib.lines import Line2D
colors = ["tab:blue", "tab:orange", "tab:red"]
fig, axes = plt.subplots(3, 3, figsize=(6, 6), sharex=True, sharey=True)
for i, col1 in enumerate(stbin.columns):
    for j, col2 in enumerate(stbin.columns):
        if i == 2:
            axes[i, j].set_xlabel(col2)
            axes[i,j].set_xticks([0,1,2],["[0, 0.1)", "[0.1, 0.9)", "[0.9, 1]"], rotation=45)
        if j == 0:
            axes[i, j].set_ylabel(col1+"\n fraction")
        if i!=j:
            tb = stbin[[stbin.columns[i], stbin.columns[j]]].value_counts().unstack().fillna(0).astype(int)
            tb  = tb/tb.sum()
            tb.index = ["[0, 0.1)", "[0.1, 0.9)", "[0.9, 1]"]
            tb.columns = ["[0, 0.1)", "[0.1, 0.9)", "[0.9, 1]"]
            tb.T.plot(kind='bar', stacked=True, ax=axes[i,j], legend=False, color=colors)
legend_categ = tb.index
legend_elements = [Line2D([0], [0], color=colors[i], lw=5, label=legend_categ[i]) for i in range(len(legend_categ))]
# Add a legend with custom legend entries
axes[0,0].legend(handles=legend_elements)
#axes[0,0].set_xlim([zmin, zmax])
axes[0,0].set_ylim([zmin, zmax])
axes[2,0].set_xticks([0,1,2],["[0, 0.1)", "[0.1, 0.9)", "[0.9, 1]"], rotation=45)
axes[2,1].set_xticks([0,1,2],["[0, 0.1)", "[0.1, 0.9)", "[0.9, 1]"], rotation=45) #somehow needed to re-do..
axes[2,2].set_xticks([0,1,2],["[0, 0.1)", "[0.1, 0.9)", "[0.9, 1]"], rotation=45)
plt.tight_layout()
plt.savefig("/Users/qingbowang/Downloads/vscoloc_corr_bar.png", dpi=500)
plt.savefig("/Users/qingbowang/Downloads/vscoloc_corr_bar.pdf", dpi=500)
plt.clf()