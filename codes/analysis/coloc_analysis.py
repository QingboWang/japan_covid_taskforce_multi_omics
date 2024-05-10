import pandas as pd

allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/multi_combined_pip0001_hg19annot.txt", sep="\t")
wb = pd.read_csv("~/Desktop/taskforce_n1102/n1300/wbfrac_classification.tsv.gz", sep="\t", compression="gzip", index_col=0)
wb.index.name = "Gene"
allj["Gene"] = allj.gene_id.str.split("\\.").str[0]
allj.set_index("Gene", inplace=True)
allj = allj.join(wb, how="left") #after setting index properly
allj.reset_index(inplace=True)
#GTEx 49 tissues PIPs:
gtex = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/GTEx_49tissues_pip0001.tsv.gz", sep="\t", compression="gzip")
gtex["Gene"] = gtex.gene.str.split("\\.").str[0]
gtex["variant_id_hg38"] = gtex.variant_hg38.str.replace("_b38","").str.replace("_",":")
allj.set_index(["variant_id_hg38", "Gene"], inplace=True)
gtex = gtex.pivot_table(index = ["variant_id_hg38","Gene"], columns = "tissue", values = "gtex_min_pip")
gtex.columns = "gtex_" + gtex.columns
allj = allj.join(gtex, how="left")
allj.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip0001_allqtl_gtexannot_forcoloc.tsv.gz", sep="\t", compression='gzip')#Go till 0.001
#read, in local:
allj = pd.read_csv("~/Downloads/jctf_pip0001_allqtl_gtexannot_forcoloc.tsv.gz", sep="\t")

#Plot for 49 tissues:
#x = 09, 01, 001 threshold
#y = eQTL colocalization with 1. JCTF / 2. GTEx WB 3. GTEx liver 4. GTEx Spleen 5. GTEx any other
#denom = num. genes with both measurement = 2,211 subset genes
genes = pd.read_csv("~/Downloads/n998_mrna.protein_regressedout.bed", sep='\t', usecols=[0,1,2,3]).gene_id.str.split("\\.").str[0]
allj.set_index("Gene", inplace=True) #For later convenience
allj['clpp_jctf'] = allj.e_min_pip * allj.p_min_pip #1. JCTF
allj['clpp_jctf_max'] = np.maximum(allj.e_pip_sus, allj.e_pip_fm) * np.maximum(allj.p_pip_sus, allj.p_pip_fm) #also the max
allj['clpp_jctf_susie'] = allj.e_pip_sus * allj.p_pip_sus
allj['clpp_jctf_fm'] = allj.e_pip_fm * allj.p_pip_fm #also the max
allj['clpp_jctf_mean'] = allj[["e_pip_sus", "e_pip_fm"]].mean(axis=1) * allj[["p_pip_sus", "p_pip_fm"]].mean(axis=1)
allj['clpp_gtex_wb'] = allj.gtex_Whole_Blood * allj.p_min_pip #2. GTEx WB
allj['clpp_gtex_liver'] = allj.gtex_Liver * allj.p_min_pip #3. GTEx Liver
allj['clpp_gtex_spleen'] = allj.gtex_Spleen * allj.p_min_pip#4. GTEx Spleen
allj["pip_gtex_max_others"] = allj.loc[:,(allj.columns.str.startswith("gtex"))&(allj.columns!="gtex_Whole_Blood")&(allj.columns!="gtex_Liver")&(allj.columns!="gtex_Spleen")].fillna(0).max(axis=1)
allj['clpp_gtex_max_others'] = allj.pip_gtex_max_others * allj.p_min_pip#5. GTEx all others
cols_keep = allj.columns[allj.columns.str.contains("clpp")]
clpp_st = allj.groupby("Gene")[cols_keep].max().fillna(0)
del clpp_st["clpp"] #duplicated column, identical as clpp_jctf
#And first, also this to impute the 0s (i.e. genes with PIP<0.01 in all measures)
rows_to_add = np.setdiff1d(genes, clpp_st.index) #filling the genes not included in allj (i.e. everything <0.01)
for row in rows_to_add:
    clpp_st.loc[row,:] = 0
clpp_st = clpp_st.loc[clpp_st.index.intersection(genes),:] #filtering to original 2211 genes in intersection
clpp_st.to_csv("/Users/qingbowang/Downloads/pqtl_clpp_stats.tsv", sep="\t")
#plot:
clpp_st = pd.read_csv("/Users/qingbowang/Downloads/pqtl_clpp_stats.tsv", sep="\t", index_col=0)
#1. x = threshold, y = either clpp min, mean or max, etc
st = {}
stmean = {}
stsem = {}
#thresholds = [0.001]+list(np.arange(0.01,1.01, 0.01))
thresholds = np.logspace(-3,0,31) #log linear
for thres in thresholds:
    st[thres] = (clpp_st>thres).sum(axis=0)
    stmean[thres] = (clpp_st > thres).mean(axis=0)
    stsem[thres] = (clpp_st > thres).sem(axis=0)
st = pd.DataFrame(st)
stmean = pd.DataFrame(stmean)
stsem = pd.DataFrame(stsem)
#Plot:
rorder = ["clpp_jctf_max", "clpp_jctf_mean", "clpp_jctf_susie", "clpp_jctf_fm", "clpp_jctf"]
leg = ["Max(SuSiE, FINEMAP)", "Mean(SuSiE, FINEMAP)", "SuSiE", "FINEMAP", "Min(SuSiE, FINEMAP)\n(our choice)"]
colors = ["tab:green", "tab:blue", "tab:brown", "tab:purple", "tab:orange"]
x = stsem.columns
plt.figure(figsize=(6,4.2))
for i in range(5):
    plt.errorbar(x, stmean.loc[rorder[i],:], stsem.loc[rorder[i],:], color=colors[i], fmt="o", label=leg[i])
plt.xlim([0.00091,1.1])
plt.xscale("log")
plt.xticks(thresholds[::10], thresholds[::10], fontsize=14)
plt.xlabel("Colocalization Posterior Probability (CLPP) threshold", fontsize=14)
plt.ylabel("Fraction of genes with colocalization\nsignal above threshold", fontsize=14)
plt.legend(title="Method:", loc="upper right")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/clpp_method_comparison.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/clpp_method_comparison.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Next, adding different method by cumulative
clpp_st["j_only"]  = clpp_st.clpp_jctf
clpp_st["j_gwb"]  = np.maximum(clpp_st.j_only, clpp_st.clpp_gtex_wb)#jctf + gtex wb
clpp_st["j_gwb_gliv"]  = np.maximum(clpp_st.j_gwb, clpp_st.clpp_gtex_liver)
clpp_st["j_gwb_gliv_gspl"]  = np.maximum(clpp_st.j_gwb_gliv, clpp_st.clpp_gtex_spleen)
clpp_st["j_gwb_gliv_gspl_oth"]  = np.maximum(clpp_st.j_gwb_gliv_gspl, clpp_st.clpp_gtex_max_others)
clpp_st = clpp_st.loc[:,clpp_st.columns.str.startswith("j_")] #the ones we care
st = {}
stmean = {}
stsem = {}
thresholds = np.logspace(-3,0,31) #log linear
for thres in thresholds:
    st[thres] = (clpp_st>thres).sum(axis=0)
    stmean[thres] = (clpp_st > thres).mean(axis=0)
    stsem[thres] = (clpp_st > thres).sem(axis=0)
st = pd.DataFrame(st)
stmean = pd.DataFrame(stmean)
stsem = pd.DataFrame(stsem)
#Plot:
leg = ["JCTF (Whole Blood)", "+GTEx Whole Blood", "+GTEx Liver", "+GTEx Spleen", "+GTEx all others"]
import matplotlib.cm as cm
vir = cm.viridis
colors = [vir(0), vir(0.2), vir(0.4), vir(0.6), vir(0.8)]
x = stsem.columns
plt.figure(figsize=(6,4.2))
for i in range(5):
    plt.errorbar(x, stmean.iloc[i,:], stsem.iloc[i,:], color=colors[i], fmt="o", label=leg[i])
plt.xlim([0.00091,1.1])
plt.xscale("log")
plt.xticks(thresholds[::10], thresholds[::10], fontsize=14)
plt.xlabel("Colocalization Posterior Probability (CLPP) threshold", fontsize=14)
plt.ylabel("Fraction of genes with colocalization\nsignal above threshold", fontsize=14)
plt.legend(title="eQTL data:", loc="upper right")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/clpp_gtex_comparison.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/clpp_gtex_comparison.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Finally, take the coloc probability per WB expression bin:
#i.e. x = CLPP threhold, y = %genes, per gene classification.
allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/multi_combined_pip0001_hg19annot.txt", sep="\t")
wbfrac = pd.read_csv("~/Desktop/taskforce_n1102/n1300/wbfrac_classification.tsv.gz", sep="\t", compression="gzip", index_col=0)
allj.index = allj.gene_id.str.split("\\.").str[0]
allj = allj.join(wbfrac.wb_frac_bin, how="left")
allj.wb_frac_bin.fillna("isNA").value_counts() #Let's probably just remove those isNA guys
allj = allj[~allj.wb_frac_bin.isna()]
#and make max per gene
allj = allj.sort_values(by=["gene_id", "clpp"], ascending=False).groupby("gene_id").head(1)
#And take the coloc per thres:
stmeanhigh = {}
stsemhigh = {}
stmeanmid = {}
stsemmid = {}
stmeanlow = {}
stsemlow = {}
thresholds = np.logspace(-3,0,31) #log linear
dthigh = allj[allj.wb_frac_bin=="high"]
dtmid = allj[allj.wb_frac_bin=="mid"]
dtlow = allj[allj.wb_frac_bin=="low"]
for thres in thresholds:    
    stmeanhigh[thres] = (dthigh.clpp > thres).mean(axis=0)
    stsemhigh[thres] = (dthigh.clpp > thres).sem(axis=0)
    stmeanmid[thres] = (dtmid.clpp > thres).mean(axis=0)
    stsemmid[thres] = (dtmid.clpp > thres).sem(axis=0)
    stmeanlow[thres] = (dtlow.clpp > thres).mean(axis=0)
    stsemlow[thres] = (dtlow.clpp > thres).sem(axis=0)
stmeanhigh = pd.DataFrame(stmeanhigh, index=["high"])
stsemhigh = pd.DataFrame(stsemhigh, index=["high"])
stmeanmid = pd.DataFrame(stmeanmid, index=["mid"])
stsemmid = pd.DataFrame(stsemmid, index=["mid"])
stmeanlow = pd.DataFrame(stmeanlow, index=["low"])
stsemlow = pd.DataFrame(stsemlow, index=["low"])
stmean = pd.concat([stmeanhigh, stmeanmid, stmeanlow])
stsem = pd.concat([stsemhigh, stsemmid, stsemlow]) #Beatiful.

#Plot:
rorder = ["high", "mid", "low"]
leg = ["High (5%$\leq$)", "Middle (1~5%)", "Low (1%$\leq$)"]
pallete = sns.cubehelix_palette(5, start=2)
colors = [pallete[-2], pallete[-3], pallete[-5]] #前classifyで使ったやつ.
x = stsem.columns
plt.figure(figsize=(6,4.2))
for i in range(3):
    plt.errorbar(x, stmean.loc[rorder[i],:], stsem.loc[rorder[i],:], color=colors[i], fmt="o", label=leg[i])
plt.xlim([0.00091,1.1])
plt.xscale("log")
plt.xticks(thresholds[::10], thresholds[::10], fontsize=14)
plt.xlabel("Colocalization Posterior Probability (CLPP) threshold", fontsize=14)
plt.ylabel("Fraction of genes with colocalization\nsignal above threshold", fontsize=14)
plt.legend(title="%Whole Blood expression", loc="upper right")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/coloc_per_wbfrac.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/coloc_per_wbfrac.pdf', bbox_inches='tight', dpi=500)
plt.clf()



#And also the effect size concordance when colocalizing
allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eps_n998_pip001.tsv.gz", sep="\t")
allj.set_index(["variant_id_hg38", "gene_id"], inplace=True)
eadd = []
padd = []
enoadd = []
pnoadd = []
for chr in range(1,23):
    e = pd.read_csv("~/Desktop/taskforce_n1102/n1300/n998_eqtl_sumstats/eqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr),sep='\t')
    p = pd.read_csv("~/Desktop/taskforce_n1102/n1300/n998_pqtl_sumstats/pqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr),sep='\t')
    e["variant_id_hg38"] = e.variant_id.str.split("_").str[0]
    p["variant_id_hg38"] = p.variant_id.str.split("_").str[0]
    e.set_index(["variant_id_hg38", "gene_id"], inplace=True)
    p.set_index(["variant_id_hg38", "gene_id"], inplace=True)
    eno = e.loc[e.index.difference(allj.index),:]
    pno = p.loc[p.index.difference(allj.index), :]
    enoadd.append(eno)
    pnoadd.append(pno)
    e = e.loc[e.index.intersection(allj.index),:]
    p = p.loc[p.index.intersection(allj.index), :]
    eadd.append(e)
    padd.append(p)
    print ("Done chr{0}, {1}".format(chr, tm.ctime()))
eadd = pd.concat(eadd)
padd = pd.concat(padd)
enoadd = pd.concat(enoadd)
pnoadd = pd.concat(pnoadd) #for 0,0 bin
enouse = enoadd.sample(n=10000, random_state=0)
pnouse = pnoadd.loc[enouse.index,:]
#enoadd.sort_index(inplace=True)
#pnoadd.sort_index(inplace=True)
allj = allj.join(eadd[["slope", "slope_se"]], how="left").join(padd[["slope", "slope_se"]], how="left", rsuffix="_p")
#and bin the PIP
def bin_pip(pip):
    if pip<0.01: return 0
    elif pip<0.1: return 0.01
    elif pip<0.9: return 0.1
    elif pip<1.1: return 0.9
    else: return np.nan
allj = allj[(allj.e_pip_bin>=0.01)|(allj.p_pip_bin>=0.01)]
allj["e_pip_bin"] = allj.e_min_pip.apply(lambda x: bin_pip(x))
allj["p_pip_bin"] = allj.p_min_pip.apply(lambda x: bin_pip(x))
tb = allj.groupby(["e_pip_bin", "p_pip_bin"]).size().unstack()
tb
#and take the effect size correlation per PIP bin:
from scipy import stats
pips = ["[0,0.01)", "[0.01,0.1)", "[0.1,0.9)", "[0.9,1]"]
colors = ["tab:blue", "tab:green", "tab:orange", "tab:red"]
#m,M = min(allj.slope.min(), allj.slope_p.min()), max(allj.slope.max(),allj.slope_p.max()) #roughly -3 to 4
fig, axes = plt.subplots(4, 4, figsize=(7.5, 7.5), sharex=True, sharey=True)
for i in range(len(tb.index)): #i for pQTL index, y axis
    for j in range(len(tb.index)): #j for eQTL index, x axis
        axes[i, j].axvline(x=0, linestyle="--", linewidth=0.5, zorder=-2, color="black")
        axes[i, j].axhline(y=0, linestyle="--", linewidth=0.5, zorder=-2, color="black")
        #axes[i,j].plot([m*1.1,M*1.1], [m*1.1,M*1.1], linestyle="--", linewidth=0.5, zorder=-2, color="black")
        axes[i, j].plot([-3.2, 4.2], [-3.2, 4.2], linestyle="--", linewidth=0.5, zorder=-2, color="black")
        if max(i,j)==0:
            x = enouse.slope
            y = pnouse.slope
            xerr = enouse.slope_se
            yerr = pnouse.slope_se
            #axes[i,j].scatter(x, y, c="tab:grey")
            axes[i, j].errorbar(x=x, y=y, yerr=yerr, xerr=xerr, color="tab:grey", alpha=0.75, mec="black", fmt="o")
            axes[i,j].text(x=-3.1, y=4.1, va="top", ha="left", s="r={0}\nn=random 10k".format(np.round(stats.pearsonr(x,y)[0],3)))
        if max(i, j) > 0:
            pipi = tb.index[i]
            pipj = tb.index[j]
            color = colors[min(i,j)]
            df = allj[(allj.e_pip_bin==pipj)&(allj.p_pip_bin==pipi)]
            x = df.slope
            y = df.slope_p
            xerr = df.slope_se
            yerr = df.slope_se_p
            #axes[i,j].scatter(x, y, c=color)
            axes[i, j].errorbar(x=x, y=y, yerr=yerr, xerr=xerr, color=color, alpha=0.75, mec="black", fmt="o")
            axes[i,j].text(x=-2.9, y=3.9, va="top", ha="left", s="r={0}\nn={1}".format(np.round(stats.pearsonr(x,y)[0],3), len(x)))
        if j==0:
            axes[i,j].set_ylabel("PIP: {0}\npQTL effect size".format(pips[i]))
        if i==3:
            axes[i,j].set_xlabel("PIP: {0}\neQTL effect size".format(pips[j]))
axes[0,0].set_xlim([-3.2,4.2])
axes[0,0].set_ylim([-3.2,4.2])
axes[0,0].set_xticks([-2,0,2,4], [-2,0,2,4])
axes[0,0].set_yticks([-2,0,2,4], [-2,0,2,4])
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/coloc_effsize_comparison.png', bbox_inches='tight', dpi=500)
#plt.savefig('/Users/qingbowang/Downloads/coloc_effsize_comparison.pdf', bbox_inches='tight', dpi=500)
plt.clf()


