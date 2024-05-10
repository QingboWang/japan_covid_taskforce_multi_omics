import pandas as pd
import numpy as np
import pickle
from scipy import stats
import glob
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"


#first read the mpra tier data:
fns = glob.glob("/Users/qingbowang/Desktop/mpra_analysis/out_dir_b1_k562_20221022/K562_alpha_pval_20221201_*")
k562_p = []
for fn in fns:
    k562_p.append(pd.read_csv(fn, sep="\t", index_col=0))
k562_p = pd.concat(k562_p)

fns = glob.glob("/Users/qingbowang/Desktop/mpra_analysis/out_dir_b1_20220909/HEPG2_alpha_pval_20221201_*")
hepg2_p = []
for fn in fns:
    hepg2_p.append(pd.read_csv(fn, sep="\t", index_col=0))
hepg2_p = pd.concat(hepg2_p)
k562_p.sort_values(by="pval", ascending=True, inplace=True)
k562_p["idx"] = np.arange(k562_p.shape[0])+1
k562_p['bh_value_001'] = k562_p.idx/sum(~k562_p.pval.isna()) * 0.01
fdr001thres = k562_p[k562_p.pval<k562_p.bh_value_001].tail(1).pval.values[0]
k562_p['bh_value_01'] = k562_p.idx/sum(~k562_p.pval.isna()) * 0.1
fdr01thres = k562_p[k562_p.pval<k562_p.bh_value_01].tail(1).pval.values[0]
k562_p["emvar"] = "None"
k562_p.loc[k562_p.pval<=fdr01thres,"emvar"] = "Tier 2"
k562_p.loc[k562_p.pval<=fdr001thres,"emvar"] = "Tier 1"
hepg2_p.sort_values(by="pval", ascending=True, inplace=True)
hepg2_p["idx"] = np.arange(hepg2_p.shape[0])+1
hepg2_p['bh_value_001'] = hepg2_p.idx/sum(~hepg2_p.pval.isna()) * 0.01
fdr001thres = hepg2_p[hepg2_p.pval<hepg2_p.bh_value_001].tail(1).pval.values[0]
hepg2_p['bh_value_01'] = hepg2_p.idx/sum(~hepg2_p.pval.isna()) * 0.1
fdr01thres = hepg2_p[hepg2_p.pval<hepg2_p.bh_value_01].tail(1).pval.values[0]
hepg2_p["emvar"] = "None"
hepg2_p.loc[hepg2_p.pval<=fdr01thres,"emvar"] = "Tier 2"
hepg2_p.loc[hepg2_p.pval<=fdr001thres,"emvar"] = "Tier 1"


#then, the list of n1019 variants (updated):
v1019 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/unique_variants_n1019.tsv", sep="\t")
v1019.index = "b1_24K_"+v1019.iloc[:,-1].str.split("_").str[0].str.replace(":","_")
v1019["exists_in_jctf"] = True
k562_p = k562_p.join(v1019.exists_in_jctf, how="left")
k562_p["exists_in_jctf"] = k562_p.exists_in_jctf.fillna(False)
hepg2_p = hepg2_p.join(v1019.exists_in_jctf, how="left")
hepg2_p["exists_in_jctf"] = hepg2_p.exists_in_jctf.fillna(False)

#and the pips:
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=' ')
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep='\t')
e_fm["variant_id_hg38"] = e_fm.rsid.str.split("_").str[0]
e_fm["gene_id"] = e_fm.rsid.str.split("_").str[-1]
e_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
e_sus["variant_id_hg38"] = e_sus.rsid.str.split("_").str[0]
e_sus["gene_id"] = e_sus.rsid.str.split("_").str[-1]
e_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
ej = e_fm.join(e_sus.pip, how="inner")
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej.columns = ["rsid", "pip_fm", "pip_sus"]
ej["min_pip"] = np.minimum(ej.pip_fm, ej.pip_sus)
ej = ej.groupby(ej.index.get_level_values(0))["min_pip"].max()
ej.index = "b1_24K_"+ej.index.str.replace(":","_")

k562_p = k562_p.join(ej, how="left")
k562_p["min_pip"] = k562_p.min_pip.fillna(0)#0 for either not tested, or below the threshold
hepg2_p = hepg2_p.join(ej, how="left")
hepg2_p["min_pip"] = hepg2_p.min_pip.fillna(0)#0 for either not tested, or below the threshold

def bin_pip(x): #returning lower bound
    if x<0.1: return (0)
    #elif x<0.1: return (0.01)
    elif x<0.5: return (0.1)
    elif x<0.9: return (0.5)
    elif x<1: return (0.9)
    else: return (1)

k562_p["pip_bin"] = k562_p.min_pip.apply(lambda x: bin_pip(x))
st_k562 = k562_p[k562_p.exists_in_jctf].groupby(["emvar", "pip_bin"]).size().unstack().fillna(0).astype(int)
st_k562 = st_k562.iloc[[0,2,1],:] #sort manually
stnorm_k562 = (st_k562/st_k562.sum())

hepg2_p["pip_bin"] = hepg2_p.min_pip.apply(lambda x: bin_pip(x))
st_hepg2 = hepg2_p[hepg2_p.exists_in_jctf].groupby(["emvar", "pip_bin"]).size().unstack().fillna(0).astype(int)
st_hepg2 = st_hepg2.iloc[[0,2,1],:] #sort manually
stnorm_hepg2 = (st_hepg2/st_hepg2.sum())

#either tissue:
#classification: (Tier 1,Tier 1), (Tier 1, Tier 2) in either order, (Tier 2, Tier 2), (Tier 2, None) in either order, and None
totalp = k562_p.join(hepg2_p.emvar, how="left", rsuffix="_hepg2")
totalp["tier"] = "None"
totalp.loc[((totalp.emvar=="Tier 2") | (totalp.emvar_hepg2=="Tier 2")) & ((totalp.emvar!="None") | (totalp.emvar_hepg2!="None")),"tier"] = "Tier 2"
totalp.loc[(totalp.emvar=="Tier 1") | (totalp.emvar_hepg2=="Tier 1"),"tier"] = "Tier 1"
#save:
totalp.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/mpra_results_summary_for_jctf_pqtl_manuscript.tsv", sep="\t")
st = totalp.groupby(["tier", "pip_bin"]).size().unstack().fillna(0).astype(int)
st = st.iloc[[0,2,1],:] #sort manually
stnorm = (st/st.sum())
#Also the error bars:
sterr_k562 = np.sqrt(stnorm_k562*(1-stnorm_k562)/st_k562.sum(axis=0))
sterr_hepg2 = np.sqrt(stnorm_hepg2*(1-stnorm_hepg2)/st_hepg2.sum(axis=0))
sterr = np.sqrt(stnorm*(1-stnorm)/st.sum(axis=0))
#plot
import matplotlib.cm as cm
vir = cm.viridis
c5 = ["tab:blue", "tab:olive", "tab:orange", "tab:red", "magenta"]
labs = ["<0.1","[0.1,0.5)", "[0.5,0.9)", "[0.9,1)", "1"]
plt.figure(figsize=(10,3.5))
plt.bar(np.arange(5), stnorm_k562.iloc[2,:]*100, color=vir(0.9), width=0.6) #tier 1
plt.bar(np.arange(5), stnorm_k562.iloc[1,:]*100, bottom=stnorm_k562.iloc[2,:]*100, color=vir(0.5), width=0.6)
plt.bar(np.arange(5)+6, stnorm_hepg2.iloc[2,:]*100, color=vir(0.9), width=0.6) #tier 1
plt.bar(np.arange(5)+6, stnorm_hepg2.iloc[1,:]*100, bottom=stnorm_hepg2.iloc[2,:]*100, color=vir(0.5), width=0.6)
plt.bar(np.arange(5)+12, stnorm.iloc[2,:]*100, color=vir(0.9), width=0.6, label="Tier 1") #tier 1
plt.bar(np.arange(5)+12, stnorm.iloc[1,:]*100, bottom=stnorm.iloc[2,:]*100, color=vir(0.5), width=0.6, label="Tier 2")
plt.errorbar(np.arange(5), stnorm_k562.iloc[2,:]*100, sterr_k562.iloc[2,:]*100, fmt="None", color="black")
plt.errorbar(np.arange(5), (stnorm_k562.iloc[1,:]+stnorm_k562.iloc[2,:])*100, sterr_k562.iloc[1,:]*100, fmt="None", color="black")
plt.errorbar(np.arange(5)+6, stnorm_hepg2.iloc[2,:]*100, sterr_hepg2.iloc[2,:]*100, fmt="None", color="black")
plt.errorbar(np.arange(5)+6, (stnorm_hepg2.iloc[1,:]+stnorm_hepg2.iloc[2,:])*100, sterr_hepg2.iloc[1,:]*100, fmt="None", color="black")
plt.errorbar(np.arange(5)+12, stnorm.iloc[2,:]*100, sterr.iloc[2,:]*100, fmt="None", color="black")
plt.errorbar(np.arange(5)+12, (stnorm.iloc[1,:]+stnorm.iloc[2,:])*100, sterr.iloc[1,:]*100, fmt="None", color="black")
plt.xticks(list(np.arange(5))+list(np.arange(5)+6)+list(np.arange(5)+12), labs+labs+labs, rotation=45)
plt.xlabel("PIP bin")
plt.ylabel("%Emvar")
plt.axvline(x=5, linestyle="--", linewidth=0.5, color="tab:gray")
plt.axvline(x=11, linestyle="--", linewidth=0.5, color="tab:gray")
plt.text(-0.5,15.5, "Cell type: K562", va="top")
plt.text(5.5,15.5, "Cell type: HepG2", va="top")
plt.text(11.5,15.5, "Cell type:\nK562 or HepG2", va="top")
plt.legend(title="MPRA emvar\nclassification", loc="upper left", bbox_to_anchor=(1, 1.1))
plt.xlim([-0.9,16.9])
for i in range(5):
    plt.gca().get_xticklabels()[i].set_color(c5[i])
    plt.gca().get_xticklabels()[i+5].set_color(c5[i])
    plt.gca().get_xticklabels()[i+10].set_color(c5[i])
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/mpra_mainfig.png", dpi=500)
plt.savefig("/Users/qingbowang/Desktop/plots/mpra_mainfig.pdf", dpi=500)
plt.clf()

#eff size:
vlist = k562_p.index.intersection(hepg2_p.index).str.replace("b1_24K_","") #list of the variant that we care
ej["variant_hg38"] = ej.rsid.str.split("_").str[0]
ej.index = ej["variant_hg38"].str.replace(":","_")
ejneed = ej.loc[ej.index.intersection(vlist),:]
ejneed["variant_id"] = ejneed.rsid.str.split("_ENSG").str[0]
ejneed["gene_id"] = ejneed.rsid.str.split("_").str[-1]
ejneed.set_index(["variant_id", "gene_id"], inplace=True) #PIPs for the variants we care, when exists

dfkeeps = []
for chr in list(range(1,23))+["X"]:
    df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr),sep='\t')
    df.index = df.variant_id.str.split("_").str[0]
    df = df.loc[df.index.intersection(vlist.str.replace("_", ":")),:]
    df.set_index(["variant_id", "gene_id"], inplace=True)
    df = df.join(ejneed.min_pip, how="left")
    df["min_pip"] = df.min_pip.fillna(0)
    #Then, filter to best PIP per variant.
    df = df.sort_values(by=["min_pip", "pval_nominal"], ascending=[False, True]).groupby("variant_id").head(1)
    dfkeeps.append(df)
    print ("done {0}, {1}".format(chr, tm.ctime()))
    print("kept {0} variant-genes".format(df.shape[0]))
dfkeeps = pd.concat(dfkeeps)
dfkeeps.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/mpra_variants_gtex_effsize.tsv", sep="\t")
dfkeeps = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/mpra_variants_gtex_effsize.tsv", sep="\t")
dfkeeps.index = "b1_24K_" + dfkeeps.variant_id.str.split("_").str[0].str.replace(":","_")
k562_p = k562_p.join(dfkeeps[["slope","slope_se"]], how="left")
hepg2_p = hepg2_p.join(dfkeeps[["slope","slope_se"]], how="left")
#and save:
k562_p.to_csv("/Users/qingbowang/Desktop/mpra_analysis/out_dir_b1_k562_20221022/K562_emvar_pvals_and_jctfpip_updated.tsv.gz", sep="\t")
hepg2_p.to_csv("/Users/qingbowang/Desktop/mpra_analysis/out_dir_b1_20220909/HEPG2_emvar_pvals_and_jctfpip_updated.tsv.gz", sep="\t")

k562_p = pd.read_csv("/Users/qingbowang/Desktop/mpra_analysis/out_dir_b1_k562_20221022/K562_emvar_pvals_and_jctfpip_updated.tsv.gz", sep="\t", index_col=0)
hepg2_p = pd.read_csv("/Users/qingbowang/Desktop/mpra_analysis/out_dir_b1_20220909/HEPG2_emvar_pvals_and_jctfpip_updated.tsv.gz", sep="\t", index_col=0)


#and plot:
def bin_pip_loose(x): #returning lower bound
    if x<0.1: return (0)
    elif x < 0.99: return (0.9)
    else: return (1)
k562_p["pip_bin_loose"] = k562_p.min_pip.apply(lambda x: bin_pip_loose(x))
k562_p["same_eff_dir"] = (k562_p.alpha_diff*k562_p.slope) > 0
dir_k562 = k562_p[(~(k562_p.alpha_diff*k562_p.slope).isna())].groupby(["emvar", "pip_bin_loose"]).same_eff_dir.agg(["mean", "sem", "count"])
hepg2_p["pip_bin_loose"] = hepg2_p.min_pip.apply(lambda x: bin_pip_loose(x))
hepg2_p["same_eff_dir"] = (hepg2_p.alpha_diff*hepg2_p.slope) > 0
dir_hepg2 = hepg2_p[(~(hepg2_p.alpha_diff*hepg2_p.slope).isna())].groupby(["emvar", "pip_bin_loose"]).same_eff_dir.agg(["mean", "sem", "count"])

#plot everything:
stnorm_k562 = dir_k562["mean"].unstack()
sterr_k562 = dir_k562["sem"].unstack()
stnorm_hepg2 = dir_hepg2["mean"].unstack()
sterr_hepg2 = dir_hepg2["sem"].unstack()

import matplotlib.cm as cm
vir = cm.viridis
c3 = ["tab:blue", "tab:orange", "magenta"]
labs = ["<0.1","[0.1,0.99)", "[0.99,1]"]
fig, ax = plt.subplots(1,2,figsize=(8,3.5), sharex=True, sharey=True)
ax[0].errorbar([-0.1,1-0.1,2-0.1], stnorm_k562.iloc[0,:]*100,
                      sterr_k562.iloc[0,:]*100, fmt="o", color="tab:gray", label="None")
ax[0].errorbar([0,1,2], stnorm_k562.iloc[2,:]*100,
                      sterr_k562.iloc[2,:]*100, fmt="o", color=vir(0.5), label="Tier 2")
ax[0].errorbar([+0.1,1+0.1,2+0.1], stnorm_k562.iloc[1,:]*100,
                      sterr_k562.iloc[1,:]*100, fmt="o", color=vir(0.9), label="Tier 1")
ax[1].errorbar([-0.1,1-0.1,2-0.1], stnorm_hepg2.iloc[0,:]*100,
                      sterr_hepg2.iloc[0,:]*100, fmt="o", color="tab:gray", label="None")
ax[1].errorbar([0,1,2], stnorm_hepg2.iloc[2,:]*100,
                      sterr_hepg2.iloc[2,:]*100, fmt="o", color=vir(0.5), label="Tier 2")
ax[1].errorbar([+0.1,1+0.1,2+0.1], stnorm_hepg2.iloc[1,:]*100,
                      sterr_hepg2.iloc[1,:]*100, fmt="o", color=vir(0.9), label="Tier 1")
fig.text(0.45, 0.01, 'PIP bin', ha='center', fontsize=16)
ax[0].set_ylabel("%Same effect direction")
ax[0].text(0,99, "Cell type: K562", va="top")
ax[1].text(0,99, "Cell type: HepG2", va="top")
ax[1].legend(title="MPRA emvar\nclassification", loc="upper left", bbox_to_anchor=(1, 1.1))
ax[0].set_ylim([46,104])
ax[0].yaxis.set_ticks([50,60,70,80,90,100])
for i in range(2):
    ax[i].set_xticks(list(np.arange(3)), labs, rotation=45)
    ax[i].grid(axis='y', color='0.95')
    ax[i].spines[['right', 'top']].set_visible(False)
for i in range(3):
    ax[0].xaxis.get_ticklabels()[i].set_color(c3[i])
    ax[1].xaxis.get_ticklabels()[i].set_color(c3[i])
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/mpra_effsize_mainfig.png", dpi=500)
plt.savefig("/Users/qingbowang/Desktop/plots/mpra_effsize_mainfig.pdf", dpi=500)
plt.clf()

