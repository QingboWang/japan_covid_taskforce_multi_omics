#all the sQTL related analysis comes here, as well as plots
import pandas as pd
import numpy as np
import time as tm
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
c6 = ["tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
from matplotlib.colors import LogNorm
import seaborn as sns
from scipy import stats

#read the data
sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/sqtl_redo_susie_pip0001.txt", sep="\t")
fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/sqtl_redo_finemap_pip0001.txt", sep=' ')
#let's first make it per variant
sus["variant_hg38"] = sus.iloc[:,0].str.split("_").str[0]
fm["variant_hg38"] = fm.iloc[:,0].str.split("_").str[0]
susv = sus.groupby("variant_hg38").pip.max()
fmv = fm.groupby("variant_hg38").prob.max()
#correlation between two PIPs:
#check 1: SuSiE vs FM for each.
ej = pd.concat([susv, fmv], axis=1)
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej.columns = ["pip_fm", "pip_sus"]
def bin_pip(x): #returning lower bound
    if x<0.001: return (0)
    elif x<0.01: return (0.001)
    elif x<0.1: return (0.01)
    elif x < 0.5: return (0.1)
    elif x<0.9: return (0.5)
    elif x<1: return (0.9)
    else: return (1)
ej["fm_bin"] = ej.pip_fm.apply(lambda x: bin_pip(x))
ej["sus_bin"] = ej.pip_sus.apply(lambda x: bin_pip(x))
ej["min_pip"] = np.minimum(ej.pip_fm, ej.pip_sus)
ej["min_pip_bin"] = ej.min_pip.apply(lambda x: bin_pip(x))
#table:
tb = ej.groupby(["fm_bin", "sus_bin"]).size().unstack().fillna(0).astype(int) #very good agreement confirmed (so far)

#Number compared to the previous release:
n_new = tb.iloc[-2,-1]+tb.iloc[-1,-2]+tb.iloc[-1,-1]
prev_rel = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_rna_releasedata_cis_sqtls.tsv.gz", sep='\t', compression='gzip')
n_old = len(prev_rel[((prev_rel.pip_fm>0.9) & (prev_rel.pip_susie>0.9))].variant_id_hg38.unique())
print (n_new, n_old) #counting the number of unique variant hg38
tb.iloc[0,0] = -1 #we do not look at this cell.
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/sqtl_sus_vs_fm_bins.tsv", sep='\t')

#Plot this tb:
tk = ["[0,0.001)","[0.001,0.01)","[0.01,0.1)","[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/sqtl_sus_vs_fm_bins.tsv", sep='\t', index_col=0)
tb.index = tk
tb.columns = tk
tb = tb.iloc[::-1,:] #to match the order
log_norm = LogNorm()
plt.figure(figsize=(8,8))
sns.heatmap(tb+1, annot=tb, fmt="d", square=True, linewidths=.5, norm=log_norm,
            cmap="viridis", cbar_kws={'label': 'count',
                                      "shrink": .6})
plt.yticks(rotation=20, fontsize=15)
plt.xticks(rotation=20, fontsize=15)
plt.xlabel("SuSiE PIP", fontsize=16)
plt.ylabel("FINEMAP PIP", fontsize=16)
plt.title("sQTL, per variant counts")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/sqtl_pip_fm_vs_susie_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/sqtl_pip_fm_vs_susie_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()


#VEP:
ej = ej[ej.min_pip>=0.001]
vps = []
for chr in list(range(1,23))+["X"]:
    vp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/all_n1419_vepped/n1419_vepped_chr{0}.txt".format(chr), sep='\t')
    vp["variant_id_hg38"] = vp["#Uploaded_variation"].str.split("_").str[0]
    vp = vp[vp.Consequence.str.contains("splice")]
    vps.append(vp)
    if chr==1:
        print (vp)
    print ("done chr{0}, {1}".format(chr, tm.ctime()))
vps = pd.concat(vps)
vps.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_splice.tsv.gz", sep='\t', compression='gzip')

#per variant stats:
univs = []
for chr in list(range(1,22+1))+["X"]:
    df = pd.read_csv("~/Desktop/taskforce_n1102/n1300/sqtl_sumstats/sqtl_n1019_redo.chr{0}.allpairs.mac2.intron5em8.txt.gz".format(chr),sep='\t', usecols=["variant_id"], squeeze=True)
    univs.append(pd.Series(df.unique()))
univs = pd.concat(univs)
univs.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/splice_uniquevars_tested.tsv.gz", sep='\t', compression='gzip')
N0 = len(univs) #is 4513469

#Get the enrichment
vps = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_splice.tsv.gz", sep='\t', index_col=0)
#make it per variant, any anotations concat by ,
vps["variant"] = vps.iloc[:,0].str.split("_").str[0]
vps.index = vps.variant
vps.index.names = ["variant_hg38"]
vps = vps.groupby('variant').agg({'Consequence': lambda x: ','.join(x)})

ej = ej.join(vps, how="left")
ej["Consequence"].fillna("NA", inplace=True)
#get the list of splice related annotations:
annots = pd.Series(",".join(vps.tolist()).split(",")).unique()
annots = pd.Series(annots)[pd.Series(annots).str.contains("splice")]
univs = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/splice_uniquevars_tested.tsv.gz", sep='\t', compression='gzip', index_col=0)
vexists = univs.iloc[:,0].str.split("_").str[0] #to get the n0, filter by variants existing in sQTL analysis at all
N0 = len(univs) #is 4513469
tb = {}
for annot in annots:
    ej[annot] = ej.Consequence.str.contains(annot)
    t = ej.groupby(["min_pip_bin", annot]).size().unstack().T.fillna(0).astype(int)
    oth = t.sum(axis=1)
    n0 = sum(vps.loc[np.intersect1d(vps.index, vexists),"Consequence"].str.contains(annot))
    t[0] = [N0-n0, n0]
    t[0] = t[0] - oth #What if we do not do this subtraction? Not real but just to see
    t = t.T
    t["frac"] = t[True]/t.sum(axis=1)
    t["err"] = np.sqrt(t.frac*(1-t.frac)/t.sum(axis=1))
    t.sort_index(inplace=True)
    tb[annot] = t
    print (t)
    t.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/sqtl_enr_{0}.tsv".format(annot), sep='\t')

#Number for splice donor and acceptor for example (and the enrichment and the fisher):
annot = "splice_donor_variant"
t1 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/sqtl_enr_{0}.tsv".format(annot), sep='\t', index_col=0)
x0 = t1.iloc[:-2,:2].sum()
x1 = t1.iloc[-2:,:2].sum()
stats.fisher_exact([x0, x1])
enr = (x1[1]/x1.sum()) /(x0[1]/x0.sum())
print (enr)
annot = "splice_acceptor_variant"
t2 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/sqtl_enr_{0}.tsv".format(annot), sep='\t', index_col=0)
x0 = t2.iloc[:-2,:2].sum()
x1 = t2.iloc[-2:,:2].sum()
stats.fisher_exact([x0, x1])
enr = (x1[1]/x1.sum()) /(x0[1]/x0.sum())
print (enr)

#plot:
tk = ["[0,0.001)","[0.001,0.01)","[0.01,0.1)","[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
ys = []
yerrs = []
for annot in annots:
    f0 = tb[annot][True].sum() / tb[annot].iloc[:,:2].sum().sum()
    ys.append(tb[annot].frac / f0)
    yerrs.append(tb[annot].err / f0)
ys = pd.DataFrame(ys)
ys.index = annots
yerrs = pd.DataFrame(yerrs)
yerrs.index = annots
ys.sort_values(by=0.9, inplace=True)
yerrs = yerrs.loc[ys.index, :]
ys_use = ys.iloc[2:,:]
yerrs_use = yerrs.iloc[2:,:]
x = np.arange(ys_use.shape[0])
plt.rcParams.update({'font.size': 14})
c7 = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
plt.figure(figsize=(7.5,4))
i = 0
for pip in ys_use.columns:
    y = ys_use[pip]
    yerr = yerrs_use[pip]
    plt.errorbar(x-0.3+0.1*i, y, yerr, fmt="o", label=tk[i], color=c7[i])
    i += 1
plt.xticks(x, ys_use.index.str.replace("_","\n"), rotation=0)
plt.xlabel("Annotation")
plt.ylabel("Enrichment")
plt.legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="sQTL PIP bin")
plt.yscale('log')
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/sqtl_vep_enr_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/sqtl_vep_enr_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#PIP vs any splice annotation:
ej["splice_evidence"] = ej.iloc[:,-6:].sum(axis=1)>0
tb = []
for i in np.arange(101) / 100:
    j = ej[ej.min_pip >= i].splice_evidence.value_counts()
    tb.append(j)
    print ("done {0}, {1}".format(i, tm.ctime()))
tb = pd.DataFrame(tb)
tb.index = range(101)
tb = tb.fillna(0).astype(int)
tb.loc[0,False] = N0-tb.loc[0,True] #replacing with the total num for 0
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/sqtl_vep_frac.tsv", sep='\t')

tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/sqtl_vep_frac.tsv", sep='\t', index_col=0)
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 2]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variants)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,1], label=tbnorm.columns[1], color="tab:orange")
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], bottom=tbnorm.iloc[:,1], label=tbnorm.columns[0], color="tab:gray")
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Has splice annotations in VEP:", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/sqtl_vepannot_frac_overview_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/sqtl_vepannot_frac_overview_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()


#sQTL enrichment for different classes of molQTLs:
allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eps_n998_pip001.tsv.gz", sep="\t")
allj["clpp"] = allj.e_min_pip*allj.p_min_pip
allj["classification"] = "NA"
allj.loc[allj.clpp>0.9,"classification"] = "coloc_strict"
allj.loc[(allj.clpp>0.1)&(allj.clpp<=0.9),"classification"] = "coloc_lenient"
allj.loc[(allj.clpp>0.01)&(allj.clpp<=0.1),"classification"] = "coloc_super_lenient"
allj.loc[(allj.ps_min_pip>0.9)&(allj.clpp<0.01),"classification"] = "protein"
allj.loc[(allj.es_min_pip>0.9)&(allj.clpp<0.01),"classification"] = "mrna"
allj.classification.value_counts()
allj.index = allj.variant_id_hg38
ej.index.names = allj.index.names
allj = allj.join(ej.min_pip_bin, rsuffix="_sqtl", how="left")
allj.rename(columns = {"min_pip_bin":"min_pip_bin_sqtl"}, inplace=True)
allj.min_pip_bin_sqtl.fillna(0,inplace=True)
f1 = allj[allj.classification=="coloc_strict"].min_pip_bin_sqtl.value_counts().sort_index()
f2 = allj[allj.classification=="coloc_lenient"].min_pip_bin_sqtl.value_counts().sort_index()
f3 = allj[allj.classification=="coloc_super_lenient"].min_pip_bin_sqtl.value_counts().sort_index()
f4 = allj[allj.classification=="protein"].min_pip_bin_sqtl.value_counts().sort_index()
f5 = allj[allj.classification=="mrna"].min_pip_bin_sqtl.value_counts().sort_index()
tb = pd.concat([f1,f2,f3,f4,f5], axis=1).sort_index().fillna(0).astype(int)
tb.columns = ["coloc_strict", "coloc_lenient", "coloc_super_lenient", "protein", "mrna"]
tb/tb.sum()
(tb/tb.sum()).iloc[::-1,:].cumsum()

tb = tb.iloc[::-1,:]
tb = tb.cumsum()
del tb["coloc_super_lenient"]
del tb["coloc_lenient"] #also no need
tb = tb.iloc[:,[0,2,1]]
tbnorm = tb/tb.iloc[-1,:]
tberr = np.sqrt((tbnorm*(1-tbnorm))/ tb.iloc[-1,:])
#the ns:
tb.sum(axis=0)
#and plot:
plt.rcParams.update({'font.size': 14})
c7 = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tk = ["[0,0.001)","[0.001,0.01)","[0.01,0.1)","[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
plt.figure(figsize=(3.5,5))
for i in range(tb.shape[0]):
    plt.bar(np.arange(tbnorm.shape[1]), tbnorm.iloc[6-i,:], color=c7[i], label=tk[i])
    plt.errorbar(np.arange(tbnorm.shape[1])+0.025*(3-i), tbnorm.iloc[6-i,:], yerr=tberr.iloc[6-i,:],
                 capsize=0, fmt="o", ms=0, color="black", elinewidth=0.8)
plt.legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="sQTL PIP bin", fontsize=12, title_fontsize=12)
plt.xlabel("QTL classification", fontsize=14)
plt.ylabel("Fraction", fontsize=14)
plt.xlim([-0.5,2.5])
plt.ylim([0,1])
plt.xticks(np.arange(tb.shape[1]),
           ["Both\n(n={0})".format(int(sum(allj.classification=="coloc_strict"))),
            "mRNA\nspecific\n(n={0})".format(int(sum(allj.classification=="mrna"))),
            "Protein\nspecific\n(n={0})".format(int(sum(allj.classification=="protein")))], rotation=30)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/epqtl_sqtlcoloc_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/epqtl_sqtlcoloc_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Fisher:
from scipy import stats
x0 = [tb.iloc[:-2,0].sum(),tb.iloc[-2:,0].sum()]
#coloc strict vs mRNA
x1 = [tb.iloc[:-2,-1].sum(),tb.iloc[-2:,-1].sum()]
stats.fisher_exact([x0,x1])
#coloc strict vs protein
x2 = [tb.iloc[:-2,-2].sum(),tb.iloc[-2:,-2].sum()]
stats.fisher_exact([x0,x2])

