
import pandas as pd
import numpy as np
import time as tm
#eQTL:
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=" ", index_col=0)
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep="\t", index_col=0)
ej = e_fm.join(e_sus.pip, how="inner")
ej.fillna(0, inplace=True)
ej.columns = ["pip_fm", "pip_sus"]
ej["pip_min"] = np.minimum(ej.pip_fm, ej.pip_sus)
ej["variant_id_hg38"] = ej.index.str.split("_").str[0]
ej["gene_id"] = ej.index.str.split("_").str[-1]
ej.set_index(["variant_id_hg38", "gene_id"], inplace=True)
#pQTL:
p_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_finemap_pip0001.txt", sep=' ')
p_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_susie_pip0001.txt", sep='\t')
p_fm["variant_id_hg38"] = p_fm.rsid.str.split("_").str[0]
p_fm["gene_id"] = p_fm.rsid.str.split("_").str[-1]
p_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
p_sus["variant_id_hg38"] = p_sus.rsid.str.split("_").str[0]
p_sus["gene_id"] = p_sus.rsid.str.split("_").str[-1]
p_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
pj = p_fm.join(p_sus.pip, how="inner")
pj.fillna(0, inplace=True) #anything lower than 0001 is 0
pj.columns = ["rsid", "pip_fm", "pip_sus"]
#Add MAF info
padd = []
eadd = []
for chr in range(1,23):
    p = pd.read_csv("~/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t')
    e = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t')
    p["variant_id_hg38"] = p.variant_id.str.split("_").str[0]
    e["variant_id_hg38"] = e.variant_id.str.split("_").str[0]
    p.set_index(["variant_id_hg38", "gene_id"], inplace=True)
    e.set_index(["variant_id_hg38", "gene_id"], inplace=True)
    p = p.loc[p.index.intersection(pj.index),["tss_distance", "ma_count", "maf"]]
    e = e.loc[e.index.intersection(ej.index), ["tss_distance", "ma_count", "maf"]]
    padd.append(p)
    eadd.append(e)
    print ("{0} p and {1} e added".format(p.shape[0], e.shape[0]))
    print ("Done chr{0}, {1}".format(chr, tm.ctime()))
padd = pd.concat(padd)
eadd = pd.concat(eadd)
padd.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip0001_pqtl_maf.tsv.gz", sep="\t", compression='gzip')
eadd.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip0001_eqtl_maf.tsv.gz", sep="\t", compression='gzip')

#eQTL:
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=" ", index_col=0)
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep="\t", index_col=0)
ej = e_fm.join(e_sus.pip, how="inner")
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej.columns = ["pip_fm", "pip_sus"]
ej["pip_min"] = np.minimum(ej.pip_fm, ej.pip_sus)
ej["variant_id_hg38"] = ej.index.str.split("_").str[0]
ej["gene_id"] = ej.index.str.split("_").str[-1]
ej.set_index(["variant_id_hg38", "gene_id"], inplace=True)
eadd = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip0001_eqtl_maf.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
ej = ej.join(eadd, how="left")
ej = ej[~ej.maf.isna()]
#pQTL:
p_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_finemap_pip0001.txt", sep=' ')
p_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_susie_pip0001.txt", sep='\t')
p_fm["variant_id_hg38"] = p_fm.rsid.str.split("_").str[0]
p_fm["gene_id"] = p_fm.rsid.str.split("_").str[-1]
p_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
p_sus["variant_id_hg38"] = p_sus.rsid.str.split("_").str[0]
p_sus["gene_id"] = p_sus.rsid.str.split("_").str[-1]
p_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
pj = p_fm.join(p_sus.pip, how="inner")
pj.fillna(0, inplace=True) #anything lower than 0001 is 0
pj.columns = ["rsid", "pip_fm", "pip_sus"]
#pj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/vepannot_pip0001_aminoacids_uniprotannot_updated.tsv.gz", sep='\t')
#pj.set_index(["variant_id_hg38", "gene_id"], inplace=True) #These are only missense, filtered.
padd = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip0001_pqtl_maf.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
pj = pj.join(padd, how="left")
#and add missense annotation:
mis = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/vepannot_pip0001_aminoacids_uniprotannot_updated.tsv.gz", sep='\t')
mis.set_index(["variant_id_hg38", "gene_id"], inplace=True) #These are only missense, filtered.
pj = pj.join(mis.missense_variant, how="left")
pj = pj[~pj.maf.isna()]
pj["missense_variant"].fillna(False, inplace=True)
#And do the enrichment analysis:
def bin_pip(pip):
    if pip<0.01: return 0
    elif pip<0.1: return 0.01
    elif pip < 0.5: return 0.1
    elif pip < 0.7: return 0.5
    elif pip<0.9: return 0.7
    elif pip<1.1: return 0.9
    else: return np.nan
#eQTL:
ej["tss_proximal"] = abs(ej.tss_distance)<5*10**2
ej.tss_proximal.value_counts()
ej["eqtl_pip_bin"] = ej.pip_min.apply(lambda x: bin_pip(x))
ej.groupby(["eqtl_pip_bin", "tss_proximal"]).size().unstack().fillna(0).astype(int)
#Stratify by MAC
ej["mac_bin"] = 1000 #1000 as the largest
ej.loc[ej.ma_count<=300,"mac_bin"] = 300
ej.loc[ej.ma_count<=100,"mac_bin"] = 100
ej.loc[ej.ma_count<=30,"mac_bin"] = 30
ej.loc[ej.ma_count<=10,"mac_bin"] = 10
ej.mac_bin.value_counts().sort_index()
y = ej.groupby(["eqtl_pip_bin", "mac_bin"]).tss_proximal.agg("mean").unstack().fillna(0)
yerr = ej.groupby(["eqtl_pip_bin", "mac_bin"]).tss_proximal.agg("sem").unstack().fillna(0)
yline = ej[ej.ma_count>=100].groupby(["eqtl_pip_bin"]).tss_proximal.agg("mean")
#Plot:
c6 = ["tab:grey", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red"]
labels = ["[0,0.01)", "[0.01,0.1)", "[0.1,0.5)", "[0.5,0.7)", "[0.7,0.9)", "[0.9,1]"]
x = np.arange(yerr.shape[1])
plt.figure(figsize=(6.5,3.5))
for i in range(yerr.shape[0]):
    plt.errorbar(x-0.3+i*0.1, y.iloc[i,:], yerr.iloc[i,:], color=c6[i], fmt="o", label=labels[i])
    if i==yerr.shape[0]-1: #Just for legend, expected to be overlayed
        plt.axhline(y=yline.values[i], color="black", linestyle="--", linewidth=0.5, zorder=-2, label="Common\nvariants'\naverage")
    plt.axhline(y=yline.values[i], color=c6[i], linestyle="--", linewidth=0.5, zorder=-2)
plt.xlabel("Minor Allele Count (MAC) of eQTL variants")
plt.ylabel("Fraction of TSS-proximal variants")
plt.xticks(x, ["(2,10]", "(10,30]", "(30,100]", "(100,300]", "300<"], rotation=30)
plt.legend(bbox_to_anchor=(1,1.05))
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_eqtl_funcenr.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_eqtl_funcenr.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#pQTL, missense var enrichment:
pj["pip_min"] = np.minimum(pj.pip_sus, pj.pip_fm)
pj["pqtl_pip_bin"] = pj.pip_min.apply(lambda x: bin_pip(x))
pj.groupby(["pqtl_pip_bin", "missense_variant"]).size().unstack().fillna(0).astype(int)
pj["mac_bin"] = 1000 #1000 as the largest
pj.loc[pj.ma_count<=300,"mac_bin"] = 300
pj.loc[pj.ma_count<=100,"mac_bin"] = 100
pj.loc[pj.ma_count<=30,"mac_bin"] = 30
pj.loc[pj.ma_count<=10,"mac_bin"] = 10
pj.mac_bin.value_counts().sort_index()
y = pj.groupby(["pqtl_pip_bin", "mac_bin"]).missense_variant.agg("mean").unstack().fillna(0)
yerr = pj.groupby(["pqtl_pip_bin", "mac_bin"]).missense_variant.agg("sem").unstack().fillna(0)
yline = pj[pj.ma_count>=100].groupby(["pqtl_pip_bin"]).missense_variant.agg("mean")
#Plot:
c6 = ["tab:grey", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red"]
labels = ["[0,0.01)", "[0.01,0.1)", "[0.1,0.5)", "[0.5,0.7)", "[0.7,0.9)", "[0.9,1]"]
x = np.arange(yerr.shape[1])
plt.figure(figsize=(6.5,3.5))
for i in range(yerr.shape[0]):
    plt.errorbar(x-0.3+i*0.1, y.iloc[i,:], yerr.iloc[i,:], color=c6[i], fmt="o", label=labels[i])
    if i==yerr.shape[0]-1: #Just for legend, expected to be overlayed
        plt.axhline(y=yline.values[i], color="black", linestyle="--", linewidth=0.5, zorder=-2, label="Common\nvariants'\naverage")
    plt.axhline(y=yline.values[i], color=c6[i], linestyle="--", linewidth=0.5, zorder=-2)
plt.xlabel("Minor Allele Count (MAC) of pQTL variants")
plt.ylabel("Fraction of missense variants")
plt.xticks(x, ["(2,10]", "(10,30]", "(30,100]", "(100,300]", "300<"], rotation=30)
plt.legend(bbox_to_anchor=(1,1.05))
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_pqtl_funcenr.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_pqtl_funcenr.pdf', bbox_inches='tight', dpi=500)
plt.clf()

pj["tss_proximal"] = abs(pj.tss_distance)<5*10**2
y = pj.groupby(["pqtl_pip_bin", "mac_bin"]).tss_proximal.agg("mean").unstack().fillna(0)
yerr = pj.groupby(["pqtl_pip_bin", "mac_bin"]).tss_proximal.agg("sem").unstack().fillna(0)
yline = pj[pj.ma_count>=100].groupby(["pqtl_pip_bin"]).tss_proximal.agg("mean")
#Plot:
c6 = ["tab:grey", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red"]
labels = ["[0,0.01)", "[0.01,0.1)", "[0.1,0.5)", "[0.5,0.7)", "[0.7,0.9)", "[0.9,1]"]
x = np.arange(yerr.shape[1])
plt.figure(figsize=(6.5,3.5))
for i in range(yerr.shape[0]):
    plt.errorbar(x-0.3+i*0.1, y.iloc[i,:], yerr.iloc[i,:], color=c6[i], fmt="o", label=labels[i])
    if i==yerr.shape[0]-1: #Just for legend, expected to be overlayed
        plt.axhline(y=yline.values[i], color="black", linestyle="--", linewidth=0.5, zorder=-2, label="Common\nvariants'\naverage")
    plt.axhline(y=yline.values[i], color=c6[i], linestyle="--", linewidth=0.5, zorder=-2)
plt.xlabel("Minor Allele Count (MAC) of pQTL variants")
plt.ylabel("Fraction of TSS-proximal variants")
plt.xticks(x, ["(2,10]", "(10,30]", "(30,100]", "(100,300]", "300<"], rotation=30)
plt.legend(bbox_to_anchor=(1,1.05))
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_pqtl_funcenr_tss.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_pqtl_funcenr_tss.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Numbers:
enum = ej.groupby(["eqtl_pip_bin", "mac_bin"]).size().unstack().fillna(0).astype(int)
enum.columns = ["(2,10]", "(10,30]", "(30,100]", "(100,300]", "300<"]
esum = enum.sum(axis=1)
efrac = (enum.T/esum).T

pnum = pj.groupby(["pqtl_pip_bin", "mac_bin"]).size().unstack().fillna(0).astype(int)
pnum.columns = ["(2,10]", "(10,30]", "(30,100]", "(100,300]", "300<"]
psum = pnum.sum(axis=1)
pfrac = (pnum.T/psum).T

#Plot:
import matplotlib.cm as cm
vir = cm.viridis
colors = [vir(0), vir(0.2), vir(0.4), vir(0.6), vir(0.8)]
plt.figure(figsize=(5,3.5))
pfrac.plot.bar(stacked=True, color=colors)
plt.xticks(np.arange(6), ["[0,0.01)\nn={0}".format(psum.values[0]),
                          "[0.01,0.1)\nn={0}".format(psum.values[1]),
                          "[0.1,0.5)\nn={0}".format(psum.values[2]),
                          "[0.5,0.7)\nn={0}".format(psum.values[3]),
                          "[0.7,0.9)\nn={0}".format(psum.values[4]),
                          "[0.9,1]\nn={0}".format(psum.values[5])], rotation=30)
handles, labels = plt.gca().get_legend_handles_labels()
handles.reverse()
labels.reverse()
plt.legend(handles, labels, bbox_to_anchor=(1, 1.02), title="MAC:")
plt.xlabel("pQTL PIP bin")
plt.ylabel("MAC distribution")
plt.text(x=4.8, y=0.12, s="n=\n{0}".format(pnum.iloc[-1,0]), color="#ffffff")
plt.ylim([0,1])
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_pqtl_numbers.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_pqtl_numbers.pdf', bbox_inches='tight', dpi=500)
plt.clf()

plt.figure(figsize=(5,3.5))
pfrac.plot.bar(stacked=True, color=colors)
plt.xticks(np.arange(6), ["[0,0.01)\nn={0}".format(esum.values[0]),
                          "[0.01,0.1)\nn={0}".format(esum.values[1]),
                          "[0.1,0.5)\nn={0}".format(esum.values[2]),
                          "[0.5,0.7)\nn={0}".format(esum.values[3]),
                          "[0.7,0.9)\nn={0}".format(esum.values[4]),
                          "[0.9,1]\nn={0}".format(esum.values[5])], rotation=30)
handles, labels = plt.gca().get_legend_handles_labels()
handles.reverse()
labels.reverse()
plt.legend(handles, labels, bbox_to_anchor=(1, 1.02), title="MAC:")
plt.xlabel("eQTL PIP bin")
plt.ylabel("MAC distribution")
plt.text(x=4.8, y=0.12, s="n=\n{0}".format(enum.iloc[-1,0]), color="#ffffff")
plt.ylim([0,1])
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_eqtl_numbers.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/rare_eqtl_numbers.pdf', bbox_inches='tight', dpi=500)
plt.clf()
