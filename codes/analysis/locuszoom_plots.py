#all the manhattan plots we use in the main figures:
import pandas as pd
import numpy as np
import time as tm
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
c4 = ["tab:gray", "tab:blue", "tab:orange", "tab:red"]
c6 = ["tab:blue", "mediumseagreen", "tab:green", "tab:olive", "tab:orange", "tab:red"]
c8 = ["tab:gray", "tab:brown", "tab:blue", "mediumseagreen", "tab:green", "tab:olive", "tab:orange", "tab:red"]
from matplotlib.colors import LogNorm
import seaborn as sns
import sys
import os
from os.path import exists
import io
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},sep='\t').rename(columns={'#CHROM': 'CHROM'})

#for annotating the genes
gns0 = pd.read_csv("~/Desktop/resources/gene_positions_for_locuszoom.txt", sep="\t")
gns0 = gns0.drop_duplicates(subset=["Gene stable ID"])


#NDST1
gene_id = "ENSG00000070614.15"
variant_id = "chr5:150488763:G:A"
gene_name = "NDST1"
tss = 150508131
tes = 150558211
cds_start = 150521255
cds_end = 150553332


#for gene name annotation later:
gns = gns0[gns0["Chromosome/scaffold name"]=="5"]
gns = gns[abs(gns['Transcription start site (TSS)']-150488763)<10**4]#l, rに合わせる.
gns = gns[~gns["Gene name"].isna()]
gns.sort_values(by="Gene start (bp)", inplace=True, ascending=True)
gns.reset_index(inplace=True)

#create LD:
"""memo: do this in shell
l=$(gunzip -c /Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.bed.gz | grep ENSG00000070614.15 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)
fn=/Users/qingbowang/Desktop/taskforce_n1102/n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_"$chr"_reIDed.vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/n1019_vcf_"$gn"_updated.tsv
#cut the file
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix -h $fn "$chr":$(($tss-1000000))-$(($tss+1000000)) > $fn_tmp
"""
n998 = pd.read_csv("~/Desktop/taskforce_n1102/n1300/pqtl_n998.combined_covariates_idadded.txt", sep="\t", index_col=0,nrows=2).columns
df_dir = "/Users/qingbowang/Desktop/tmp/n1019_vcf_{0}_updated.tsv".format(gene_id)
z_pqtl_dir="/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gene_id)
df = read_vcf(df_dir) #from eps_complex_trait.py
z = pd.read_csv(z_pqtl_dir, sep=' ')
df.index = df.ID #since it is already well re-IDed
#set the order same (which should be, by default)
df = df.loc[z.rsid.str.split("_ENSG").str[0],:]
#and turn it to dosage matrix;
X0 = df.iloc[:,9:].applymap(lambda x: float(x.split(":")[1])).T
X = X0.loc[n998, z.rsid.str.split("_ENSG").str[0]]
print ("done getting dosage mat for gene{0}: {1}".format(gene_id, tm.ctime()))
print ("Starting creating X^tX/N, {0}".format(tm.ctime()))
X = X.apply(lambda l: (l - l.mean())/l.std(), axis=0)
n = X.shape[0]#number of samples
R = X.T @ X / n
print ("done calculating LD mat for gene{0}: {1}".format(gene_id, tm.ctime()))
print (R.shape)
print (R.iloc[:5,:5])
#and write as .ld:
R = np.array(R)
print ("starting writing, {0}".format(tm.ctime()))
np.savetxt("/Users/qingbowang/Desktop/tmp/ld_raw_n998_{0}.ld".format(gene_id), R, delimiter=' ', fmt='%4f')
print ("done writing, {0}".format(tm.ctime()))


#lead var +-0.0015*10**8 = 15*10**4
l = 150488763 - 15*10**4
r = 150488763 + 15*10**4
#1. manhattan plot:
z_pqtl_dir="/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gene_id)
z = pd.read_csv(z_pqtl_dir, sep=' ')
z.index = z.rsid.str.split("_").str[0]
ld = pd.read_csv("/Users/qingbowang/Desktop/tmp/ld_raw_n998_{0}.ld".format(gene_id), sep=" ", header=None)
pos = z.position.astype(int)
ld.index = z.index
ld.columns = ld.index
ld = ld.loc[(l<pos)&(pos<r),variant_id]**2
#pqtl marginal:
chr = 5
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/n998_pqtl_sumstats/pqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
df = pd.read_csv(input_file, sep="\t")
df = df[df.gene_id==gene_id]
df.index = df.variant_id.str.split("_").str[0]
df = df.loc[ld.index, :]
#eqtl marginal:
chr = 5
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/n998_eqtl_sumstats/eqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
dfe = pd.read_csv(input_file, sep="\t")
dfe = dfe[dfe.gene_id==gene_id]
dfe.index = dfe.variant_id.str.split("_").str[0]
dfe = dfe.loc[ld.index, :]
#e and p pips:
allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/multi_combined_pip0001_hg19annot.txt", sep="\t")
pips = allj[allj.gene_id==gene_id]
pips["position"] = pips.variant_id_hg38.str.split(":").str[1].astype(int)
pips = pips.loc[(l<pips.position)&(pips.position<r),:]

#plot:
df["pos"] = df.index.str.split(":").str[1].astype(int)
dfe["pos"] = dfe.index.str.split(":").str[1].astype(int)
fig, ax = plt.subplots(6, 1, sharex=True, sharey=False, figsize=(8, 6), gridspec_kw={'height_ratios': [2.5,1,0.04, 2.5,1,0.3]})
ax[0].scatter(df[ld>=0].pos, -1*np.log10(df[ld>=0].pval_nominal), color="darkblue", edgecolor="black", label="[0,0.2]")
ax[0].scatter(df[ld>0.2].pos, -1*np.log10(df[ld>0.2].pval_nominal), color="lightblue", edgecolor="black", label="(0.2,0.4]")
ax[0].scatter(df[ld>0.4].pos, -1*np.log10(df[ld>0.4].pval_nominal), color="tab:green", edgecolor="black", label="(0.4,0.6]")
ax[0].scatter(df[ld>0.6].pos, -1*np.log10(df[ld>0.6].pval_nominal), color="tab:orange", edgecolor="black", label="(0.6,0.8]")
ax[0].scatter(df[ld>0.8].pos, -1*np.log10(df[ld>0.8].pval_nominal), color="tab:red", edgecolor="black", label="(0.8,0.1]")
ax[0].scatter(df[df.index==variant_id].pos, -1*np.log10(df[df.index==variant_id].pval_nominal), color="#ffffffff")
ax[0].scatter(df[df.index==variant_id].pos, -1*np.log10(df[df.index==variant_id].pval_nominal), color="tab:purple", edgecolor="black", marker="D")
ax[0].legend(title="$r^2$", fontsize=10, title_fontsize=11, loc="upper left")
ax[1].scatter(pips.position, pips.p_min_pip, color="tab:olive", edgecolor="black")
ax[1].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].p_min_pip, color="#ffffffff")
ax[1].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].p_min_pip, color="tab:olive", edgecolor="black", marker="D")
ax[3].scatter(dfe[ld>=0].pos, -1*np.log10(dfe[ld>=0].pval_nominal), color="darkblue", edgecolor="black", label="[0,0.2]")
ax[3].scatter(dfe[ld>0.2].pos, -1*np.log10(dfe[ld>0.2].pval_nominal), color="lightblue", edgecolor="black", label="(0.2,0.4]")
ax[3].scatter(dfe[ld>0.4].pos, -1*np.log10(dfe[ld>0.4].pval_nominal), color="tab:green", edgecolor="black", label="(0.4,0.6]")
ax[3].scatter(dfe[ld>0.6].pos, -1*np.log10(dfe[ld>0.6].pval_nominal), color="tab:orange", edgecolor="black", label="(0.6,0.8]")
ax[3].scatter(dfe[ld>0.8].pos, -1*np.log10(dfe[ld>0.8].pval_nominal), color="tab:red", edgecolor="black", label="(0.8,0.1]")
ax[3].legend(title="$r^2$", fontsize=10, title_fontsize=11, loc="upper left")
ax[3].text(dfe[dfe.index==variant_id].pos.values[0]+10**3.5, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal.values[0])*0.95, "rs12519827 (upstream variant)")
ax[3].scatter(dfe[dfe.index==variant_id].pos, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal), color="#ffffffff")
ax[3].scatter(dfe[dfe.index==variant_id].pos, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal), color="tab:purple", edgecolor="black", marker="D")
ax[4].scatter(pips.position, pips.e_min_pip, color="tab:blue", edgecolor="black")
ax[4].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].e_min_pip, color="#ffffffff")
ax[4].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].e_min_pip, color="tab:blue", edgecolor="black", marker="D")
ax[1].set_ylim([-0.11, 1.11])
ax[4].set_ylim([-0.11, 1.11])
ax[0].set_title(gene_name)
ax[0].set_ylabel("pQTL -log$_{10}$(p)")
ax[1].set_ylabel("pQTL\nPIP")
ax[3].set_ylabel("eQTL  -log$_{10}$(p)")
ax[4].set_ylabel("eQTL\nPIP")
ax[2].set_visible(False)
for i in [0,1,3,4]:
    ax[i].axvline(x=150488763, color="black", linestyle="--", linewidth=0.5, zorder=-1)
for i in [1,4]:
    ax[i].spines.right.set_visible(False)
    ax[i].spines.top.set_visible(False)
ax[0].set_ylim([-10,240])
#ax[0].set_ylim([-1,24])
ax[3].set_ylim([-10,240])
#ax[5] for the annotation of the gene of interest
#ax[5].plot([gns[gns["Gene name"] == gene_name]["Gene start (bp)"], gns[gns["Gene name"] == gene_name]["Gene end (bp)"]],[1, 1],color = "purple", linewidth = 4)  # the one of interest
ax[5].plot([tss,tes],[1, 1],color = "purple", linewidth = 2)  # the one of interest
ax[5].plot([cds_start,cds_end],[1, 1],color = "purple", linewidth = 4)  # the one of interest
ax[5].text(x=(tss+tes)/2, y=3, s=r"$NDST1\rightarrow$", color="purple", horizontalalignment='center',verticalalignment='center')
ax[5].set_ylim([0,4.5])
ax[5].set_yticks([])
ax[5].set_xlabel("Position on chr5 (bp)", fontsize=14 )
ax[5].spines.left.set_visible(False)
ax[5].spines.right.set_visible(False)
ax[5].spines.top.set_visible(False)
plt.tight_layout()
plt.subplots_adjust(hspace=.08)
plt.savefig("/Users/qingbowang/Desktop/plots/{0}_pips_pubfig_updated_test.png".format(gene_id), dpi=500)
plt.savefig("/Users/qingbowang/Desktop/plots/{0}_pips_pubfig_updated_test.pdf".format(gene_id), dpi=500)
plt.clf()

#TNFRSF11A
gene_id = "ENSG00000141655.17"
variant_id = "chr18:62387624:A:C"
gene_name = "TNFRSF11A"
#lead var +-0.005*10**7 = 5*10**4
l = 62387624 - 5*10**4
r = 62387624 + 5*10**4
chr = 18
#create LD:
"""memo: do thi in shell
l=$(gunzip -c /Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.bed.gz | grep ENSG00000141655.17 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)
fn=/Users/qingbowang/Desktop/taskforce_n1102/n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_"$chr"_reIDed.vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/n1419_vcf_"$gn"_updated.tsv
#cut the file
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix -h $fn "$chr":$(($tss-1000000))-$(($tss+1000000)) > $fn_tmp
"""
n998 = pd.read_csv("~/Desktop/taskforce_n1102/n1300/pqtl_n998.combined_covariates_idadded.txt", sep="\t", index_col=0,nrows=2).columns
df_dir = "/Users/qingbowang/Desktop/tmp/n1419_vcf_{0}_updated.tsv".format(gene_id)
z_pqtl_dir="/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gene_id)
df = read_vcf(df_dir) #from eps_complex_trait.py
z = pd.read_csv(z_pqtl_dir, sep=' ')
df.index = df.ID #since it is already well re-IDed
#set the order same (which should be, by default)
df = df.loc[z.rsid.str.split("_ENSG").str[0],:]
#and turn it to dosage matrix;
X0 = df.iloc[:,9:].applymap(lambda x: float(x.split(":")[1])).T
X = X0.loc[n998, z.rsid.str.split("_ENSG").str[0]]
print ("done getting dosage mat for gene{0}: {1}".format(gene_id, tm.ctime()))
print ("Starting creating X^tX/N, {0}".format(tm.ctime()))
X = X.apply(lambda l: (l - l.mean())/l.std(), axis=0)
n = X.shape[0]#number of samples
R = X.T @ X / n
print ("done calculating LD mat for gene{0}: {1}".format(gene_id, tm.ctime()))
print (R.shape)
print (R.iloc[:5,:5])
#and write as .ld:
R = np.array(R)
print ("starting writing, {0}".format(tm.ctime()))
np.savetxt("/Users/qingbowang/Desktop/tmp/ld_raw_n998_{0}.ld".format(gene_id), R, delimiter=' ', fmt='%4f')
print ("done writing, {0}".format(tm.ctime()))

#prep. manhattan plot:
z_pqtl_dir="/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gene_id)
z = pd.read_csv(z_pqtl_dir, sep=' ')
z.index = z.rsid.str.split("_").str[0]
ld = pd.read_csv("/Users/qingbowang/Desktop/tmp/ld_raw_n998_{0}.ld".format(gene_id), sep=" ", header=None)
pos = z.position.astype(int)
ld.index = z.index
ld.columns = ld.index
ld = ld.loc[(l<pos)&(pos<r),variant_id]**2
#pqtl marginal:
chr = 18
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/n998_pqtl_sumstats/pqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
df = pd.read_csv(input_file, sep="\t")
df = df[df.gene_id==gene_id]
df.index = df.variant_id.str.split("_").str[0]
df = df.loc[ld.index, :]
#eqtl marginal:
chr = 18
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/n998_eqtl_sumstats/eqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
dfe = pd.read_csv(input_file, sep="\t")
dfe = dfe[dfe.gene_id==gene_id]
dfe.index = dfe.variant_id.str.split("_").str[0]
dfe = dfe.loc[ld.index, :]
#e and p pips:
allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/multi_combined_pip0001_hg19annot.txt", sep="\t")
pips = allj[allj.gene_id==gene_id]
pips["position"] = pips.variant_id_hg38.str.split(":").str[1].astype(int)
pips = pips.loc[(l<pips.position)&(pips.position<r),:]
#bbj pips: (hg19 pos 18:60054857:A:C)
import glob
f = glob.glob("/Users/qingbowang/Desktop/resources/bbj_fm_79traits_hum0197/hum0197*/*FINEMAP*")
for fn in f:
    if fn.split("BBJ.")[1].split(".Kanai")[0]=="ALP":
        dfbbj_fm = pd.read_csv(fn, sep="\t")
f = glob.glob("/Users/qingbowang/Desktop/resources/bbj_fm_79traits_hum0197/hum0197*/*SuSiE*")
for fn in f:
    if fn.split("BBJ.")[1].split(".Kanai")[0] == "ALP":
        dfbbj_sus = pd.read_csv(fn, sep="\t")
dfbbj_fm = dfbbj_fm[(60054857-5*10**4 < dfbbj_fm.position)&(dfbbj_fm.position<60054857+5*10**4)&(dfbbj_fm.chromosome==18)]
dfbbj_sus = dfbbj_sus[(60054857-5*10**4 < dfbbj_sus.position)&(dfbbj_sus.position<60054857+5*10**4)&(dfbbj_sus.chromosome==18)]
dfbbj_fm.index = dfbbj_fm.variant
dfbbj_sus.index = dfbbj_sus.variant
dfbbj_fm = dfbbj_fm.join(dfbbj_sus.pip, rsuffix="_susie", how="outer")
dfbbj_fm["min_pip"] = np.minimum(dfbbj_fm.pip, dfbbj_fm.pip_susie)
dfbbj_fm = dfbbj_fm[dfbbj_fm.min_pip>0.00001]

#finally, the gene location annotation:
# http://asia.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000141655;r=18:62325310-62391288;t=ENST00000586569
#tss to tes = 62325310-62391288
tss = 62325310
tes = 62391288
cds_start = 62325353
cds_end = 62385034

df["pos"] = df.index.str.split(":").str[1].astype(int)
dfe["pos"] = dfe.index.str.split(":").str[1].astype(int)
fig, ax = plt.subplots(8, 1, sharex=True, sharey=False, figsize=(8, 7.8), gridspec_kw={'height_ratios': [2.5,1,0.04, 2.5,1,0.04, 1, 0.3]})
ax[0].scatter(df[ld>=0].pos, -1*np.log10(df[ld>=0].pval_nominal), color="darkblue", edgecolor="black", label="[0,0.2]")
ax[0].scatter(df[ld>0.2].pos, -1*np.log10(df[ld>0.2].pval_nominal), color="lightblue", edgecolor="black", label="(0.2,0.4]")
ax[0].scatter(df[ld>0.4].pos, -1*np.log10(df[ld>0.4].pval_nominal), color="tab:green", edgecolor="black", label="(0.4,0.6]")
ax[0].scatter(df[ld>0.6].pos, -1*np.log10(df[ld>0.6].pval_nominal), color="tab:orange", edgecolor="black", label="(0.6,0.8]")
ax[0].scatter(df[ld>0.8].pos, -1*np.log10(df[ld>0.8].pval_nominal), color="tab:red", edgecolor="black", label="(0.8,0.1]")
ax[0].scatter(df[df.index==variant_id].pos, -1*np.log10(df[df.index==variant_id].pval_nominal), color="#ffffffff")
ax[0].scatter(df[df.index==variant_id].pos, -1*np.log10(df[df.index==variant_id].pval_nominal), color="tab:purple", edgecolor="black", marker="D")
ax[0].legend(title="$r^2$", fontsize=10, title_fontsize=11, loc="upper right")
ax[0].text(df[df.index==variant_id].pos.values[0]+10**3.2, -1*np.log10(df[df.index==variant_id].pval_nominal.values[0])*0.85, "rs884205\n(3'UTR variant)")
ax[1].scatter(pips.position, pips.p_min_pip, color="tab:olive", edgecolor="black")
ax[1].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].p_min_pip, color="#ffffffff")
ax[1].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].p_min_pip, color="tab:olive", edgecolor="black", marker="D")
ax[3].scatter(dfe[ld>=0].pos, -1*np.log10(dfe[ld>=0].pval_nominal), color="darkblue", edgecolor="black", label="[0,0.2]")
ax[3].scatter(dfe[ld>0.2].pos, -1*np.log10(dfe[ld>0.2].pval_nominal), color="lightblue", edgecolor="black", label="(0.2,0.4]")
ax[3].scatter(dfe[ld>0.4].pos, -1*np.log10(dfe[ld>0.4].pval_nominal), color="tab:green", edgecolor="black", label="(0.4,0.6]")
ax[3].scatter(dfe[ld>0.6].pos, -1*np.log10(dfe[ld>0.6].pval_nominal), color="tab:orange", edgecolor="black", label="(0.6,0.8]")
ax[3].scatter(dfe[ld>0.8].pos, -1*np.log10(dfe[ld>0.8].pval_nominal), color="tab:red", edgecolor="black", label="(0.8,0.1]")
ax[3].legend(title="$r^2$", fontsize=10, title_fontsize=11, loc="upper right")
ax[3].scatter(dfe[dfe.index==variant_id].pos, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal), color="#ffffffff")
ax[3].scatter(dfe[dfe.index==variant_id].pos, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal), color="tab:purple", edgecolor="black", marker="D")
ax[4].scatter(pips.position, pips.e_min_pip, color="tab:blue", edgecolor="black")
ax[4].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].e_min_pip, color="#ffffffff")
ax[4].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].e_min_pip, color="tab:blue", edgecolor="black", marker="D")
ax[6].scatter(dfbbj_fm.position+(62387624-60054857), dfbbj_fm.min_pip, color="magenta", edgecolor="black")
ax[6].scatter(dfbbj_fm[dfbbj_fm.min_pip==1].position+(62387624-60054857), dfbbj_fm[dfbbj_fm.min_pip==1].min_pip, color="#ffffffff")
ax[6].scatter(dfbbj_fm[dfbbj_fm.min_pip==1].position+(62387624-60054857), dfbbj_fm[dfbbj_fm.min_pip==1].min_pip, color="magenta", edgecolor="black", marker="D")
#ax[6].axvspan(cds_end, tes, alpha=0.1, color='tab:gray', label="3'UTR")
#ax[6].legend(fontsize=12, loc="upper right")
ax[1].set_ylim([-0.11, 1.11])
ax[4].set_ylim([-0.11, 1.11])
ax[6].set_ylim([-0.11, 1.11])
ax[0].set_title(gene_name)
ax[0].set_ylabel("pQTL\n-log$_{10}$(p)")
ax[1].set_ylabel("pQTL\nPIP")
ax[3].set_ylabel("eQTL\n-log$_{10}$(p)")
ax[4].set_ylabel("eQTL\nPIP")
ax[6].set_ylabel("ALP PIP\n(BBJ)")
ax[2].set_visible(False)
ax[5].set_visible(False)
for i in [0,1,3,4,6]:
    #ax[i].axvspan(cds_end, tes, alpha=0.1, color='tab:gray')#no need
    ax[i].axvline(x=62387624, color="black", linestyle="--", linewidth=0.5, zorder=-1)
for i in [1,4]:
    ax[i].spines.right.set_visible(False)
    ax[i].spines.top.set_visible(False)
#for gene name annotation - wip
ax[7].plot([tss,tes],[1, 1],color = "purple", linewidth = 1.5)  # the one of interest
ax[7].plot([cds_start,cds_end],[1, 1],color = "purple", linewidth = 4)  # the one of interest
ax[7].text(x=(tss+tes)/2, y=3, s=r"$TNFRSF11A\rightarrow$", color="purple", horizontalalignment='center',verticalalignment='center')
ax[7].set_ylim([0,4.5])
ax[7].set_yticks([])
ax[7].spines.left.set_visible(False)
ax[7].spines.right.set_visible(False)
ax[7].spines.top.set_visible(False)
ax[7].set_xlabel("Position on chr18 (bp)")
ax[0].set_ylim([-0.5,14])
ax[3].set_ylim([-0.5,14])
ax[0].set_xlim(l-10**3,r+10**3)
plt.tight_layout()
plt.subplots_adjust(hspace=.08)
plt.savefig("/Users/qingbowang/Desktop/plots/{0}_pips_pubfig_updated.png".format(gene_id), dpi=500)
plt.savefig("/Users/qingbowang/Desktop/plots/{0}_pips_pubfig_updated.pdf".format(gene_id), dpi=500)
plt.clf()


#APOE
gene_id = "ENSG00000130203.10"
variant_id = "chr19:44908684:T:C"
gene_name = "APOE"
tss = 44905796
tes = 44909393
cds_start = 44906625
cds_end = 44909250

#create LD:
"""memo: do thi in shell
l=$(gunzip -c /Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.bed.gz | grep ENSG00000130203.10 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)
fn=/Users/qingbowang/Desktop/taskforce_n1102/n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_"$chr"_reIDed.vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/n1419_vcf_"$gn"_updated.tsv
#cut the file
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix -h $fn "$chr":$(($tss-1000000))-$(($tss+1000000)) > $fn_tmp
"""
n998 = pd.read_csv("~/Desktop/taskforce_n1102/n1300/pqtl_n998.combined_covariates_idadded.txt", sep="\t", index_col=0,nrows=2).columns
df_dir = "/Users/qingbowang/Desktop/tmp/n1419_vcf_{0}_updated.tsv".format(gene_id)
z_pqtl_dir="/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gene_id)
df = read_vcf(df_dir) #from eps_complex_trait.py
z = pd.read_csv(z_pqtl_dir, sep=' ')
df.index = df.ID #since it is already well re-IDed
#set the order same (which should be, by default)
df = df.loc[z.rsid.str.split("_ENSG").str[0],:]
#and turn it to dosage matrix;
X0 = df.iloc[:,9:].applymap(lambda x: float(x.split(":")[1])).T
X = X0.loc[n998, z.rsid.str.split("_ENSG").str[0]]
print ("done getting dosage mat for gene{0}: {1}".format(gene_id, tm.ctime()))
print ("Starting creating X^tX/N, {0}".format(tm.ctime()))
X = X.apply(lambda l: (l - l.mean())/l.std(), axis=0)
n = X.shape[0]#number of samples
R = X.T @ X / n
print ("done calculating LD mat for gene{0}: {1}".format(gene_id, tm.ctime()))
print (R.shape)
print (R.iloc[:5,:5])
#and write as .ld:
R = np.array(R)
print ("starting writing, {0}".format(tm.ctime()))
np.savetxt("/Users/qingbowang/Desktop/tmp/ld_raw_n998_{0}.ld".format(gene_id), R, delimiter=' ', fmt='%4f')
print ("done writing, {0}".format(tm.ctime()))

#lead var +-0.005*10**7 = 5*10**4
l = 44908684 - 5*10**4
r = 44908684 + 5*10**4
#1. manhattan plot:
z_pqtl_dir="/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gene_id)
z = pd.read_csv(z_pqtl_dir, sep=' ')
z.index = z.rsid.str.split("_").str[0]
ld = pd.read_csv("/Users/qingbowang/Desktop/tmp/ld_raw_n998_{0}.ld".format(gene_id), sep=" ", header=None)
pos = z.position.astype(int)
ld.index = z.index
ld.columns = ld.index
ld = ld.loc[(l<pos)&(pos<r),variant_id]**2
#pqtl marginal:
chr = 19
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/n998_pqtl_sumstats/pqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
df = pd.read_csv(input_file, sep="\t")
df = df[df.gene_id==gene_id]
df.index = df.variant_id.str.split("_").str[0]
df = df.loc[ld.index, :]
#eqtl marginal:
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/n998_eqtl_sumstats/eqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
dfe = pd.read_csv(input_file, sep="\t")
dfe = dfe[dfe.gene_id==gene_id]
dfe.index = dfe.variant_id.str.split("_").str[0]
dfe = dfe.loc[ld.index, :]
#e and p pips:
allj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/multi_combined_pip0001_hg19annot.txt", sep="\t")
pips = allj[allj.gene_id==gene_id]
pips["position"] = pips.variant_id_hg38.str.split(":").str[1].astype(int)
pips = pips.loc[(l<pips.position)&(pips.position<r),:]
#bbj pips: (hg19 pos 19-45411941-T-C)
import glob
f = glob.glob("/Users/qingbowang/Desktop/resources/bbj_fm_79traits_hum0197/hum0197*/*FINEMAP*")
for fn in f:
    if fn.split("BBJ.")[1].split(".Kanai")[0]=="TC":
        dfbbj_fm = pd.read_csv(fn, sep="\t")
f = glob.glob("/Users/qingbowang/Desktop/resources/bbj_fm_79traits_hum0197/hum0197*/*SuSiE*")
for fn in f:
    if fn.split("BBJ.")[1].split(".Kanai")[0] == "TC":
        dfbbj_sus = pd.read_csv(fn, sep="\t")
dfbbj_fm = dfbbj_fm[(45411941-5*10**4 < dfbbj_fm.position)&(dfbbj_fm.position<45411941+5*10**4)&(dfbbj_fm.chromosome==19)]
dfbbj_sus = dfbbj_sus[(45411941-5*10**4 < dfbbj_sus.position)&(dfbbj_sus.position<45411941+5*10**4)&(dfbbj_sus.chromosome==19)]
dfbbj_fm.index = dfbbj_fm.variant
dfbbj_sus.index = dfbbj_sus.variant
dfbbj_fm = dfbbj_fm.join(dfbbj_sus.pip, rsuffix="_susie", how="outer")
dfbbj_fm["min_pip"] = np.minimum(dfbbj_fm.pip, dfbbj_fm.pip_susie)
dfbbj_fm = dfbbj_fm[dfbbj_fm.min_pip>0.00001]

#finally, the gene location annotation (not needed in this case):
df["pos"] = df.index.str.split(":").str[1].astype(int)
dfe["pos"] = dfe.index.str.split(":").str[1].astype(int)
fig, ax = plt.subplots(8, 1, sharex=True, sharey=False, figsize=(8, 7.8), gridspec_kw={'height_ratios': [2.5,1,0.04, 2.5,1,0.04, 1, 0.3]})
ax[0].scatter(df[ld>=0].pos, -1*np.log10(df[ld>=0].pval_nominal), color="darkblue", edgecolor="black", label="[0,0.2]")
ax[0].scatter(df[ld>0.2].pos, -1*np.log10(df[ld>0.2].pval_nominal), color="lightblue", edgecolor="black", label="(0.2,0.4]")
ax[0].scatter(df[ld>0.4].pos, -1*np.log10(df[ld>0.4].pval_nominal), color="tab:green", edgecolor="black", label="(0.4,0.6]")
ax[0].scatter(df[ld>0.6].pos, -1*np.log10(df[ld>0.6].pval_nominal), color="tab:orange", edgecolor="black", label="(0.6,0.8]")
ax[0].scatter(df[ld>0.8].pos, -1*np.log10(df[ld>0.8].pval_nominal), color="tab:red", edgecolor="black", label="(0.8,0.1]")
ax[0].scatter(df[df.index==variant_id].pos, -1*np.log10(df[df.index==variant_id].pval_nominal), color="#ffffffff")
ax[0].scatter(df[df.index==variant_id].pos, -1*np.log10(df[df.index==variant_id].pval_nominal), color="tab:purple", edgecolor="black", marker="D")
ax[0].legend(title="$r^2$", fontsize=10, title_fontsize=11, loc="upper right")
ax[0].text(df[df.index==variant_id].pos.values[0]+10**3.2, -1*np.log10(df[df.index==variant_id].pval_nominal.values[0])*0.97, "rs429358 (missense)")
ax[1].scatter(pips.position, pips.p_min_pip, color="tab:olive", edgecolor="black")
ax[1].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].p_min_pip, color="#ffffffff")
ax[1].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].p_min_pip, color="tab:olive", edgecolor="black", marker="D")
ax[3].scatter(dfe[ld>=0].pos, -1*np.log10(dfe[ld>=0].pval_nominal), color="darkblue", edgecolor="black", label="[0,0.2]")
ax[3].scatter(dfe[ld>0.2].pos, -1*np.log10(dfe[ld>0.2].pval_nominal), color="lightblue", edgecolor="black", label="(0.2,0.4]")
ax[3].scatter(dfe[ld>0.4].pos, -1*np.log10(dfe[ld>0.4].pval_nominal), color="tab:green", edgecolor="black", label="(0.4,0.6]")
ax[3].scatter(dfe[ld>0.6].pos, -1*np.log10(dfe[ld>0.6].pval_nominal), color="tab:orange", edgecolor="black", label="(0.6,0.8]")
ax[3].scatter(dfe[ld>0.8].pos, -1*np.log10(dfe[ld>0.8].pval_nominal), color="tab:red", edgecolor="black", label="(0.8,0.1]")
ax[3].legend(title="$r^2$", fontsize=10, title_fontsize=11, loc="upper right")
ax[3].scatter(dfe[dfe.index==variant_id].pos, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal), color="#ffffffff")
ax[3].scatter(dfe[dfe.index==variant_id].pos, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal), color="tab:purple", edgecolor="black", marker="D")
ax[4].scatter(pips.position, pips.e_min_pip, color="tab:blue", edgecolor="black")
ax[4].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].e_min_pip, color="#ffffffff")
ax[4].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].e_min_pip, color="tab:blue", edgecolor="black", marker="D")
ax[6].scatter(dfbbj_fm.position+(44908684-45411941), dfbbj_fm.min_pip, color="magenta", edgecolor="black")
ax[6].scatter(dfbbj_fm[dfbbj_fm.min_pip==1].position+(44908684-45411941), dfbbj_fm[dfbbj_fm.min_pip==1].min_pip, color="#ffffffff")
ax[6].scatter(dfbbj_fm[dfbbj_fm.min_pip==1].position+(44908684-45411941), dfbbj_fm[dfbbj_fm.min_pip==1].min_pip, color="magenta", edgecolor="black", marker="D")
ax[1].set_ylim([-0.11, 1.11])
ax[4].set_ylim([-0.11, 1.11])
ax[6].set_ylim([-0.11, 1.11])
ax[0].set_title(gene_name)
ax[0].set_ylabel("pQTL\n-log$_{10}$(p)")
ax[1].set_ylabel("pQTL\nPIP")
ax[3].set_ylabel("eQTL\n-log$_{10}$(p)")
ax[4].set_ylabel("eQTL\nPIP")
ax[6].set_ylabel("CT PIP\n(BBJ)")
ax[2].set_visible(False)
ax[5].set_visible(False)
for i in [0,1,3,4,6]:
    ax[i].axvline(x=44908684, color="black", linestyle="--", linewidth=0.5, zorder=-1)
for i in [1,4]:
    ax[i].spines.right.set_visible(False)
    ax[i].spines.top.set_visible(False)
ax[0].set_ylim([-0.5,70])
ax[3].set_ylim([-0.5,70])
ax[0].set_xlim(l-10**3,r+10**3)
ax[7].plot([tss,tes],[1, 1],color = "purple", linewidth = 2)  # the one of interest
ax[7].plot([cds_start,cds_end],[1, 1],color = "purple", linewidth = 4)  # the one of interest
ax[7].text(x=(tss+tes)/2, y=3, s=r"$APOE\rightarrow$", color="purple", horizontalalignment='center',verticalalignment='center')
ax[7].set_ylim([0,4.5])
ax[7].set_yticks([])
ax[7].spines.left.set_visible(False)
ax[7].spines.right.set_visible(False)
ax[7].spines.top.set_visible(False)
ax[7].set_xlabel("Position on chr19 (bp)")
plt.tight_layout()
plt.subplots_adjust(hspace=.08)
plt.savefig("/Users/qingbowang/Desktop/plots/{0}_pips_pubfig_updated.png".format(gene_id), dpi=500)
plt.savefig("/Users/qingbowang/Desktop/plots/{0}_pips_pubfig_updated.pdf".format(gene_id), dpi=500)
plt.clf()


#EFHD1
gene_id = "ENSG00000115468.12"
variant_id = "chr2:232655544:A:T"
gene_name = "EFHD1"
tss = 232633604
tes = 232682776
cds_start = 232633705
cds_end = 232681719 #These are for ENST00000264059.8

#create LD:
"""memo: do thi in shell
l=$(gunzip -c /Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.bed.gz | grep ENSG00000115468.12 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)
fn=/Users/qingbowang/Desktop/taskforce_n1102/n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_"$chr"_reIDed.vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/n1419_vcf_"$gn"_updated.tsv
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix -h $fn "$chr":$(($tss-1000000))-$(($tss+1000000)) > $fn_tmp 
"""
n998 = pd.read_csv("~/Desktop/taskforce_n1102/n1300/pqtl_n998.combined_covariates_idadded.txt", sep="\t", index_col=0,nrows=2).columns
df_dir = "/Users/qingbowang/Desktop/tmp/n1419_vcf_{0}_updated.tsv".format(gene_id)
z_pqtl_dir="/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gene_id)
df = read_vcf(df_dir) #from eps_complex_trait.py
z = pd.read_csv(z_pqtl_dir, sep=' ')
df.index = df.ID #since it is already well re-IDed
#set the order same (which should be, by default)
df = df.loc[z.rsid.str.split("_ENSG").str[0],:]
#and turn it to dosage matrix;
X0 = df.iloc[:,9:].applymap(lambda x: float(x.split(":")[1])).T
X = X0.loc[n998, z.rsid.str.split("_ENSG").str[0]]
print ("done getting dosage mat for gene{0}: {1}".format(gene_id, tm.ctime()))
print ("Starting creating X^tX/N, {0}".format(tm.ctime()))
X = X.apply(lambda l: (l - l.mean())/l.std(), axis=0)
n = X.shape[0]#number of samples
R = X.T @ X / n
print ("done calculating LD mat for gene{0}: {1}".format(gene_id, tm.ctime()))
print (R.shape)
print (R.iloc[:5,:5])
#and write as .ld:
R = np.array(R)
print ("starting writing, {0}".format(tm.ctime()))
np.savetxt("/Users/qingbowang/Desktop/tmp/ld_raw_n998_{0}.ld".format(gene_id), R, delimiter=' ', fmt='%4f')
print ("done writing, {0}".format(tm.ctime()))

#lead var +-0.0005*10**7 = 5*10**4
l = 232655544 - 5*10**4
r = 232655544 + 5*10**4
#1. manhattan plot:
z_pqtl_dir="/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n998/{0}_p_fminput.z".format(gene_id)
z = pd.read_csv(z_pqtl_dir, sep=' ')
z.index = z.rsid.str.split("_").str[0]
ld = pd.read_csv("/Users/qingbowang/Desktop/tmp/ld_raw_n998_{0}.ld".format(gene_id), sep=" ", header=None)
pos = z.position.astype(int)
ld.index = z.index
ld.columns = ld.index
ld = ld.loc[(l<pos)&(pos<r),variant_id]**2
#pqtl marginal:
chr = 2
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/n998_pqtl_sumstats/pqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
df = pd.read_csv(input_file, sep="\t")
df = df[df.gene_id==gene_id]
df.index = df.variant_id.str.split("_").str[0]
df = df.loc[ld.index, :]
#eqtl marginal:
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/n998_eqtl_sumstats/eqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
dfe = pd.read_csv(input_file, sep="\t")
dfe = dfe[dfe.gene_id==gene_id]
dfe.index = dfe.variant_id.str.split("_").str[0]
dfe = dfe.loc[ld.index, :]
#e and p pips:
pips = allj[allj.gene_id==gene_id]
pips["position"] = pips.variant_id_hg38.str.split(":").str[1].astype(int)
pips = pips.loc[(l<pips.position)&(pips.position<r),:]
#bbj pips: (hg19 pos 2-233520254-A-T)
import glob
f = glob.glob("/Users/qingbowang/Desktop/resources/bbj_fm_79traits_hum0197/hum0197*/*FINEMAP*")
for fn in f:
    if fn.split("BBJ.")[1].split(".Kanai")[0]=="ALT":
        dfbbj_fm = pd.read_csv(fn, sep="\t")
f = glob.glob("/Users/qingbowang/Desktop/resources/bbj_fm_79traits_hum0197/hum0197*/*SuSiE*")
for fn in f:
    if fn.split("BBJ.")[1].split(".Kanai")[0] == "ALT":
        dfbbj_sus = pd.read_csv(fn, sep="\t")
dfbbj_fm = dfbbj_fm[(233520254-5*10**4 < dfbbj_fm.position)&(dfbbj_fm.position<233520254+5*10**4)&(dfbbj_fm.chromosome==2)]
dfbbj_sus = dfbbj_sus[(233520254-5*10**4 < dfbbj_sus.position)&(dfbbj_sus.position<233520254+5*10**4)&(dfbbj_sus.chromosome==2)]
dfbbj_fm.index = dfbbj_fm.variant
dfbbj_sus.index = dfbbj_sus.variant
dfbbj_fm = dfbbj_fm.join(dfbbj_sus.pip, rsuffix="_susie", how="outer")
dfbbj_fm["min_pip"] = np.minimum(dfbbj_fm.pip, dfbbj_fm.pip_susie)
#dfbbj_fm = dfbbj_fm[dfbbj_fm.min_pip>0.00001]

#GTEx liver PIP
cd = '/Users/qingbowang/Desktop/enformer_prep'
dfg = pd.read_csv("{0}/GTEx_49tissues_pip01.tsv.gz".format(cd), sep="\t")
dfg = dfg[dfg.tissue == "Liver"]
dfg["gene"] = dfg.gene.str.split("\\.").str[0]#unversion, for joining
dfg.set_index(["variant_hg38", "gene", "tissue"], inplace=True)
dfgs = dfg[dfg.method=="SUSIE"]
dfgf = dfg[dfg.method=="FINEMAP"]
dfg = pd.concat([dfgs.pip, dfgf.pip], axis=1).fillna(0) #missing = 0
dfg.columns = ["gtex_pip_susie", "gtex_pip_fm"]
dfg["gtex_min_pip"] = dfg.min(axis=1)
dfg["variant_hg38"] = dfg.index.get_level_values(0)
dfg = dfg[dfg.gtex_min_pip>0]
dfg = dfg[dfg.index.get_level_values(1)=="ENSG00000115468"]
dfg["position"] = dfg.variant_hg38.str.split("_").str[1].astype(int)

df["pos"] = df.index.str.split(":").str[1].astype(int)
dfe["pos"] = dfe.index.str.split(":").str[1].astype(int)
fig, ax = plt.subplots(10, 1, sharex=True, sharey=False, figsize=(8, 8.2), gridspec_kw={'height_ratios': [2.5,1,0.03, 2.5,1,0.03, 1, 0.03, 1, 0.3]})
ax[0].scatter(df[ld>=0].pos, -1*np.log10(df[ld>=0].pval_nominal), color="darkblue", edgecolor="black", label="[0,0.2]")
ax[0].scatter(df[ld>0.2].pos, -1*np.log10(df[ld>0.2].pval_nominal), color="lightblue", edgecolor="black", label="(0.2,0.4]")
ax[0].scatter(df[ld>0.4].pos, -1*np.log10(df[ld>0.4].pval_nominal), color="tab:green", edgecolor="black", label="(0.4,0.6]")
ax[0].scatter(df[ld>0.6].pos, -1*np.log10(df[ld>0.6].pval_nominal), color="tab:orange", edgecolor="black", label="(0.6,0.8]")
ax[0].scatter(df[ld>0.8].pos, -1*np.log10(df[ld>0.8].pval_nominal), color="tab:red", edgecolor="black", label="(0.8,0.1]")
ax[0].scatter(df[df.index==variant_id].pos, -1*np.log10(df[df.index==variant_id].pval_nominal), color="#ffffffff")
ax[0].scatter(df[df.index==variant_id].pos, -1*np.log10(df[df.index==variant_id].pval_nominal), color="tab:purple", edgecolor="black", marker="D")
ax[0].legend(title="$r^2$", fontsize=10, title_fontsize=11, loc="upper right")
ax[0].text(df[df.index==variant_id].pos.values[0]+10**3, -1*np.log10(df[df.index==variant_id].pval_nominal.values[0])*0.88, "rs13395911\n(intron variant)")
ax[1].scatter(pips.position, pips.p_min_pip, color="tab:olive", edgecolor="black")
ax[1].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].p_min_pip, color="#ffffffff")
ax[1].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].p_min_pip, color="tab:olive", edgecolor="black", marker="D")
ax[3].scatter(dfe[ld>=0].pos, -1*np.log10(dfe[ld>=0].pval_nominal), color="darkblue", edgecolor="black", label="[0,0.2]")
ax[3].scatter(dfe[ld>0.2].pos, -1*np.log10(dfe[ld>0.2].pval_nominal), color="lightblue", edgecolor="black", label="(0.2,0.4]")
ax[3].scatter(dfe[ld>0.4].pos, -1*np.log10(dfe[ld>0.4].pval_nominal), color="tab:green", edgecolor="black", label="(0.4,0.6]")
ax[3].scatter(dfe[ld>0.6].pos, -1*np.log10(dfe[ld>0.6].pval_nominal), color="tab:orange", edgecolor="black", label="(0.6,0.8]")
ax[3].scatter(dfe[ld>0.8].pos, -1*np.log10(dfe[ld>0.8].pval_nominal), color="tab:red", edgecolor="black", label="(0.8,0.1]")
ax[3].legend(title="$r^2$", fontsize=10, title_fontsize=11, loc="upper right")
ax[3].scatter(dfe[dfe.index==variant_id].pos, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal), color="#ffffffff")
ax[3].scatter(dfe[dfe.index==variant_id].pos, -1*np.log10(dfe[dfe.index==variant_id].pval_nominal), color="tab:purple", edgecolor="black", marker="D")
ax[3].text(dfe[dfe.index=="chr2:232636954:G:T"].pos.values[0]+10**3, -1*np.log10(dfe[dfe.index=="chr2:232636954:G:T"].pval_nominal.values[0])*0.88, "rs73995831\n(intron variant)")
ax[3].text(l+10**2, 42, "(JCTF, whole blood)")
ax[4].scatter(pips.position, pips.e_min_pip, color="tab:blue", edgecolor="black")
ax[4].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].e_min_pip, color="#ffffffff")
ax[4].scatter(pips[pips.variant_id_hg38==variant_id].position, pips[pips.variant_id_hg38==variant_id].e_min_pip, color="tab:blue", edgecolor="black", marker="D")
ax[4].text(l+10**2, 0.9, "(JCTF, whole blood)")
#ax[6].scatter(dfg.position, dfg.gtex_min_pip, color="brown", edgecolor="black")
#ax[6].scatter(dfg[dfg.variant_hg38==variant_id].position, dfg[dfg.variant_hg38==variant_id].gtex_min_pip, color="#ffffffff") #there's only one
#ax[6].scatter(dfg[dfg.variant_hg38==variant_id].position, dfg[dfg.variant_hg38==variant_id].gtex_min_pip, color="brown", edgecolor="black", marker="D")
ax[6].scatter(dfg.position, dfg.gtex_min_pip, color="brown", edgecolor="black", marker="D")
ax[6].text(l+10**2, 0.88, "(GTEx v8, liver)")
ax[8].scatter(dfbbj_fm.position+(232655544-233520254), dfbbj_fm.min_pip, color="magenta", edgecolor="black")
ax[8].scatter(dfbbj_fm[dfbbj_fm.min_pip>0.8].position+(232655544-233520254), dfbbj_fm[dfbbj_fm.min_pip>0.8].min_pip, color="#ffffffff")
ax[8].scatter(dfbbj_fm[dfbbj_fm.min_pip>0.8].position+(232655544-233520254), dfbbj_fm[dfbbj_fm.min_pip>0.8].min_pip, color="magenta", edgecolor="black", marker="D")
ax[1].set_ylim([-0.11, 1.11])
ax[4].set_ylim([-0.11, 1.11])
ax[6].set_ylim([-0.11, 1.11])
ax[8].set_ylim([-0.11, 1.11])
ax[0].set_title(gene_name)
ax[0].set_ylabel("pQTL\n-log$_{10}$(p)")
ax[1].set_ylabel("pQTL\nPIP")
ax[3].set_ylabel("eQTL\n-log$_{10}$(p)")
ax[4].set_ylabel("eQTL\nPIP")
ax[6].set_ylabel("eQTL PIP\n (GTEx v8)")
ax[8].set_ylabel("ALT PIP\n(BBJ)")
for i in [2,5,7]:
    ax[i].set_visible(False)
for i in [0,1,3,4,6,8]:
    ax[i].axvline(x=232655544, color="tab:gray", linestyle="--", linewidth=0.5, zorder=-1)
for i in [3,4,6,8]:
    ax[i].axvline(x=232636954, color="tab:gray", linestyle=":", linewidth=0.5, zorder=-1)
for i in [1, 4, 6, 8]:
    ax[i].spines.right.set_visible(False)
    ax[i].spines.top.set_visible(False)
ax[0].set_ylim([-0.5,47])
ax[3].set_ylim([-0.5,47])
ax[0].set_xlim(l-10**3,r+10**3)
#gene annotation:
ax[9].plot([tss,tes],[1, 1],color = "purple", linewidth = 2)  # the one of interest
ax[9].plot([cds_start,cds_end],[1, 1],color = "purple", linewidth = 4)  # the one of interest
ax[9].text(x=(tss+tes)/2, y=3, s=r"$EFHD1\rightarrow$", color="purple", horizontalalignment='center',verticalalignment='center')
ax[9].set_ylim([0,4.5])
ax[9].set_yticks([])
ax[9].spines.left.set_visible(False)
ax[9].spines.right.set_visible(False)
ax[9].spines.top.set_visible(False)
ax[9].set_xlabel("Position on chr2 (bp)")
plt.tight_layout()
plt.subplots_adjust(hspace=.08)
plt.savefig("/Users/qingbowang/Desktop/plots/{0}_pips_pubfig_updated.png".format(gene_id), dpi=500)
plt.savefig("/Users/qingbowang/Desktop/plots/{0}_pips_pubfig_updated.pdf".format(gene_id), dpi=500)
plt.clf()
