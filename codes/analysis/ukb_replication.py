#UKB replication


rd = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_beta.tsv.gz", sep="\t")
ps = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_beta.tsv.gz", sep="\t", compression="gzip")
pj = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001.tsv.gz", sep="\t", compression="gzip")
vp = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_and_rd10K_vep.tsv.gz", sep="\t")
afr = pd.read_csv("/home/qwang/ukb_pqtl_varnums/afr_allvars.tsv.gz", sep="\t", squeeze=True)
afr = afr.str.split(":").apply(lambda x: ":".join(x[:4]))
eur = pd.read_csv("/home/qwang/ukb_pqtl_varnums/eur_allvars.tsv.gz", sep="\t", squeeze=True)
eur = eur.str.split(":").apply(lambda x: ":".join(x[:4]))
eas = pd.read_csv("/home/qwang/ukb_pqtl_varnums/eas_allvars.tsv.gz", sep="\t", squeeze=True)
eas = eas.str.split(":").apply(lambda x: ":".join(x[:4]))
#Need to first turn everything into hg38? No actually turn ours down to hg19
rd["variant_hg19"] = rd.variant_id.str.split("_").str[-1]
rd.index = rd.variant_hg19
rd["is_in_ukb_afr"] = False
rd.loc[rd.index.intersection(afr),"is_in_ukb_afr"] = True
rd["is_in_ukb_eur"] = False
rd.loc[rd.index.intersection(eur),"is_in_ukb_eur"] = True
rd["is_in_ukb_eas"] = False
rd.loc[rd.index.intersection(eas),"is_in_ukb_eas"] = True
rd[['gene_id', 'variant_id', "is_in_ukb_afr", "is_in_ukb_eur", "is_in_ukb_eas"]].to_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_beta_ukbexists.tsv.gz", sep="\t", compression="gzip")
pj["variant_hg19"] = pj.rsid.str.split("_").str[1]
pj.index = pj.variant_hg19
pj["is_in_ukb_afr"] = False
pj.loc[pj.index.intersection(afr),"is_in_ukb_afr"] = True
pj["is_in_ukb_eur"] = False
pj.loc[pj.index.intersection(eur),"is_in_ukb_eur"] = True
pj["is_in_ukb_eas"] = False
pj.loc[pj.index.intersection(eas),"is_in_ukb_eas"] = True
pj[['rsid', "is_in_ukb_afr", "is_in_ukb_eur", "is_in_ukb_eas"]].to_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_ukbexists.tsv.gz", sep="\t", compression="gzip")
#And actually annotate the pvals:
#To do so, map ENSG to canonoical gene name:
ensgids = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/gencode.v30.genes.parsed.tsv.gz", sep='\t')
ensgids["ensgid_unv"] = ensgids.ensg_id.str.split("\\.").str[0]
ensgids.index =ensgids.ensgid_unv
#1. rd: Run this, where the python code is as below
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_ukb_data_for_rd.py eas"
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_ukb_data_for_rd.py eur"
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_ukb_data_for_rd.py afr"
"""get_ukb_data_for_rd.py
import glob
import time as tm
import pandas as pd
import numpy as np
import sys
pop = sys.argv[1]
ensgids = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/gencode.v30.genes.parsed.tsv.gz", sep='\t')
ensgids["ensgid_unv"] = ensgids.ensg_id.str.split("\\.").str[0]
ensgids.index =ensgids.ensgid_unv
rd = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_beta.tsv.gz", sep="\t")
rd["variant_hg19"] = rd.variant_id.str.split("_").str[-1]
rd["ensgid_unv"] = rd.gene_id.str.split("\\.").str[0]
rd.index = rd.ensgid_unv
rd = rd.join(ensgids["gene_name"], how="left")
gns_to_lookup = rd.gene_name.unique()
df0 = []
i = 0
files = pd.Series(glob.glob("/home/qwang/ukb_pqtl_{0}_p005/*".format(pop)))
for gn in gns_to_lookup:
    try:
        fn = files[files.str.contains("/{0}:".format(gn))].values[0]
        df = pd.read_csv(fn, sep="\t")
        df["variant_hg19"] = df.ID.str.split(":").apply(lambda x: ":".join(x[:4]))
        df.index = df.variant_hg19
        variants_to_lookup = rd[rd.gene_name==gn].variant_hg19.unique()
        df = df.loc[df.index.intersection(variants_to_lookup),:]
        df = df.iloc[:,3:-1]
        df["gene_name"] = gn
        df0.append(df)
    except:
        print ("{0} does not exist".format(gn))
    print ("Done {0}, {1}".format(i, tm.ctime()))
    i += 1
df0 = pd.concat(df0)
df0.index.name = "variant_hg19"
df0.to_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_ukb_{0}.tsv.gz".format(pop), sep="\t")
"""

#2. ps: Run this, where the python code is as below
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_ukb_data_for_ps.py eas"
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_ukb_data_for_ps.py eur"
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_ukb_data_for_ps.py afr"
"""get_ukb_data_for_ps.py
import glob
import time as tm
import pandas as pd
import numpy as np
import sys
pop = sys.argv[1]
ensgids = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/gencode.v30.genes.parsed.tsv.gz", sep='\t')
ensgids["ensgid_unv"] = ensgids.ensg_id.str.split("\\.").str[0]
ensgids.index =ensgids.ensgid_unv
ps = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_beta.tsv.gz", sep="\t", compression="gzip")
ps["variant_hg19"] = ps.variant_id.str.split("_").str[-1]
ps["ensgid_unv"] = ps.gene_id.str.split("\\.").str[0]
ps.index = ps.ensgid_unv
ps = ps.join(ensgids["gene_name"], how="left")
gns_to_lookup = ps.gene_name.unique()
df0 = []
i = 0
files = pd.Series(glob.glob("/home/qwang/ukb_pqtl_{0}_p005/*".format(pop)))
for gn in gns_to_lookup:
    try:
        fn = files[files.str.contains("/{0}:".format(gn))].values[0]
        df = pd.read_csv(fn, sep="\t")
        df["variant_hg19"] = df.ID.str.split(":").apply(lambda x: ":".join(x[:4]))
        df.index = df.variant_hg19
        variants_to_lookup = ps[ps.gene_name==gn].variant_hg19.unique()
        df = df.loc[df.index.intersection(variants_to_lookup),:]
        df = df.iloc[:,3:-1]
        df["gene_name"] = gn
        df0.append(df)
    except:
        print ("{0} does not exist".format(gn))
    print ("Done {0}, {1}".format(i, tm.ctime()))
    i += 1
df0 = pd.concat(df0)
df0.index.name = "variant_hg19"
df0.to_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_ukb_{0}.tsv.gz".format(pop), sep="\t")
"""

#Also get the lead variants, by running the code below
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_lead_var_ukb.py eas"
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_lead_var_ukb.py eur"
#qsub -l s_vmem=16G -b y "python3 ~/codes/get_lead_var_ukb.py afr"
"""get_lead_var_ukb.py
import pandas as pd
import time as tm
import glob
import sys
pop = sys.argv[1]
df0 = []
i = 0
gl_eur = pd.Series(glob.glob("/home/qwang/ukb_pqtl_{0}_p005/*".format(pop)))
for fn in gl_eur:
    df = pd.read_csv(fn, sep="\t")
    max_idx = df.LOG10P.idxmax()
    lead = df.iloc[max_idx,:]
    gn = fn.split("/")[4].split(":")[0]
    lead["gene"] = gn
    df0.append(lead)
    print ("Done {0}, {1}".format(i, tm.ctime()))
    i += 1
df0 = pd.concat(df0, axis=1).T
df0.reset_index(inplace=True)
df0.to_csv("/home/qwang/ss_for_ukb_pqtl/ukb_{0}_lead.tsv.gz".format(pop), sep="\t", compression="gzip")
"""

#Now plot the % lead var: focusing on EAS (or actually, combine the result for plots)
rd = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_beta.tsv.gz", sep="\t")
ps = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_beta.tsv.gz", sep="\t", compression="gzip")
pj = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001.tsv.gz", sep="\t", compression="gzip")
vp = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_and_rd10K_vep.tsv.gz", sep="\t")

leads = {}
for pop in ["eas", "eur", "afr"]:
    leads[pop] = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/ukb_{0}_lead.tsv.gz".format(pop), sep="\t", compression="gzip", index_col=0)
    leads[pop]["variant_hg19"] = leads[pop].ID.str.split(":").apply(lambda x: ":".join(x[:4]))
    leads[pop].set_index(["variant_hg19","gene"], inplace=True)
    leads[pop].index.names = ["variant_hg19","gene_name"]
    leads[pop]["lead_{0}".format(pop)] = True
#annotate lead var or not for each dataset (random and PIP selected):
ensgids = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/gencode.v30.genes.parsed.tsv.gz", sep='\t')
ensgids["ensgid_unv"] = ensgids.ensg_id.str.split("\\.").str[0]
ensgids.index =ensgids.ensgid_unv
rd["variant_hg19"] = rd.variant_id.str.split("_").str[-1]
rd["ensgid_unv"] = rd.gene_id.str.split("\\.").str[0]
rd.index = rd.ensgid_unv
rd = rd.join(ensgids["gene_name"], how="left")

def bin_pip(pip):
    if pip<0.001: return 0
    elif pip<0.01: return 0.001
    elif pip<0.1: return 0.01
    elif pip<0.5: return 0.1
    elif pip<0.9: return 0.5
    elif pip<1: return 0.9
    elif pip==1: return 1
    #elif pip<=1: return 0.9
    else: return np.nan
pj["p_min_pip"] = np.minimum(pj.pip_fm, pj.pip_sus)
pj["pqtl_pip_bin"] = pj.p_min_pip.apply(lambda x: bin_pip(x))
pj.pqtl_pip_bin.value_counts() #0 -> other QTL PIPが入っているだけ. can be removed
pj = pj[pj.p_min_pip>=0.001]

pj["variant_hg19"] = pj.rsid.str.split("_").str[1]
pj["ensgid_unv"] = pj.gene_id.str.split("\\.").str[0]
pj.index = pj.ensgid_unv
pj = pj.join(ensgids["gene_name"], how="left")

rd.set_index(["variant_hg19","gene_name"], inplace=True)
pj.set_index(["variant_hg19","gene_name"], inplace=True)
for pop in ["eas", "eur", "afr"]:
    rd = rd.join(leads[pop]["lead_{0}".format(pop)], how="left")
    rd["lead_{0}".format(pop)] = rd["lead_{0}".format(pop)].fillna(False)
    pj = pj.join(leads[pop]["lead_{0}".format(pop)], how="left")
    pj["lead_{0}".format(pop)] = pj["lead_{0}".format(pop)].fillna(False)
#get the stats:
#after removing chrX for fair comparison:
pj = pj[~pj.variant_id_hg38.str.startswith("chrX")]
rd = rd[~rd.variant_id_hg38.str.startswith("chrX")]
st_rd = {}
st_case = {}
for pop in ["eas", "eur", "afr"]:
    st_rd[pop] = rd["lead_{0}".format(pop)].value_counts()
    st_case[pop] = pj.groupby(["pqtl_pip_bin", "lead_{0}".format(pop)]).size().unstack().fillna(0).astype(int)
    st_case[pop].to_csv("/home/qwang/ss_for_ukb_pqtl/frac_leadvar_{0}.tsv".format(pop), sep="\t")
#also get the stats for "either"
pj["lead_any"] = pj.lead_afr | pj.lead_eur | pj.lead_eas
pj.groupby(["pqtl_pip_bin", "lead_any"]).size().unstack().fillna(0).astype(int).to_csv("/home/qwang/ss_for_ukb_pqtl/frac_leadvar_anypop.tsv", sep="\t")
#Confirmed that pd.DataFrame(st_rd) is all FALSE, no need to plot.
#Now move to local and plot these:
st = {}
frac = {}
err = {}
for pop in ["eas", "eur", "afr", "anypop"]:
    st[pop] = pd.read_csv("~/Downloads/frac_leadvar_{0}.tsv".format(pop), sep="\t", index_col=0)
    frac[pop] = (st[pop].T/st[pop].sum(axis=1)).T
    err[pop] = np.sqrt((frac[pop]*(1-frac[pop])).T/st[pop].sum(axis=1)).T
ys = pd.concat([frac["eas"]["True"], frac["eur"]["True"], frac["afr"]["True"], frac["anypop"]["True"]], axis=1)
yerrs = pd.concat([err["eas"]["True"], err["eur"]["True"], err["afr"]["True"], err["anypop"]["True"]], axis=1)
pops = ["EAS", "EUR", "AFR", "ANY"]
ys.columns = pops
yerrs.columns = pops
colors = ["tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
xlabels = [r"0.001$\leq$PIP<0.01", r"0.01$\leq$PIP<0.1", r"0.1$\leq$PIP<0.5", r"0.5$\leq$PIP<0.9", r"0.9$\leq$PIP<1","PIP=1"]
x = np.arange(len(xlabels))
plt.figure(figsize=(7,3))
for i in range(ys.shape[0]):
    plt.errorbar(x=np.arange(ys.shape[1])+i*0.1-0.25, y=ys.iloc[i,:], yerr=yerrs.iloc[i,:], color=colors[i], fmt='o', label = xlabels[i])
plt.legend(title="pQTL PIP in JCTF", bbox_to_anchor=(1.05,1))
plt.ylabel("Fraction of lead variant\nin UKB-PPP")
plt.xlabel("Population in UKB-PPP")
plt.xticks(np.arange(ys.shape[1]), pops)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/ukb_leadvar.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/ukb_leadvar.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#2. eff size concordance, as a function of PIP
#annotate rd
rd = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_beta.tsv.gz", sep="\t")
ensgids = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/gencode.v30.genes.parsed.tsv.gz", sep='\t')
ensgids["ensgid_unv"] = ensgids.ensg_id.str.split("\\.").str[0]
ensgids.index =ensgids.ensgid_unv
rd["variant_hg19"] = rd.variant_id.str.split("_").str[-1]
rd["ensgid_unv"] = rd.gene_id.str.split("\\.").str[0]
rd.index = rd.ensgid_unv
rd = rd.join(ensgids["gene_name"], how="left")
rd.set_index(["variant_hg19","gene_name"], inplace=True)
beta = {}
for pop in ["eas", "eur", "afr"]:
    beta[pop] = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_ukb_{0}.tsv.gz".format(pop), sep="\t")
    beta[pop].rename(columns = {"gene":"gene_name"}, inplace=True)
    beta[pop].set_index(["variant_hg19", "gene_name"], inplace=True)
    beta[pop].columns = beta[pop].columns + "_{0}".format(pop)
    rd = rd.join(beta[pop][["BETA"+ "_{0}".format(pop), "SE"+ "_{0}".format(pop), "LOG10P"+ "_{0}".format(pop)]], how="left")
rd.to_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_beta_ukbannot.tsv.gz", sep="\t", compression="gzip")

#annotate cases
ps = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_beta.tsv.gz", sep="\t", compression="gzip")
ps["variant_hg19"] = ps.variant_id.str.split("_").str[-1]
ps["ensgid_unv"] = ps.gene_id.str.split("\\.").str[0]
ps.index = ps.ensgid_unv
ps = ps.join(ensgids["gene_name"], how="left")
ps.set_index(["variant_hg19","gene_name"], inplace=True)
beta = {}
for pop in ["eas", "eur", "afr"]:
    beta[pop] = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_ukb_{0}.tsv.gz".format(pop), sep="\t")
    beta[pop].rename(columns = {"gene":"gene_name"}, inplace=True)
    beta[pop].set_index(["variant_hg19", "gene_name"], inplace=True)
    beta[pop].columns = beta[pop].columns + "_{0}".format(pop)
    ps = ps.join(beta[pop][["BETA"+ "_{0}".format(pop), "SE"+ "_{0}".format(pop), "LOG10P"+ "_{0}".format(pop)]], how="left")
ps.to_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_ukbannot.tsv.gz", sep="\t", compression="gzip")

#Plots:
ps = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_ukbannot.tsv.gz", sep="\t", compression="gzip")
#annotate the PIPs:
pj = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001.tsv.gz", sep="\t", compression="gzip")
ps["rsid"] = ps.variant_id + "_" + ps.gene_id
ps.set_index("rsid", inplace=True)
pj.set_index("rsid", inplace=True)
pj = pj.join(ps.iloc[:,5:], how="left")
pj.to_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_pip0001_ukbbetaannot.tsv.gz", sep="\t", compression="gzip")

#Plot eff size scatter per PIPs:
def bin_pip(pip):
    if pip<0.001: return 0
    elif pip<0.01: return 0.001
    elif pip<0.1: return 0.01
    elif pip<0.5: return 0.1
    elif pip<0.9: return 0.5
    elif pip<=1: return 0.9
    else: return np.nan
pj["p_min_pip"] = np.minimum(pj.pip_fm, pj.pip_sus)
pj["pqtl_pip_bin"] = pj.p_min_pip.apply(lambda x: bin_pip(x))
pj.pqtl_pip_bin.value_counts() #0 -> other QTL PIPが入っているだけ. can be removed
pj = pj[pj.p_min_pip>=0.001]

rd = pd.read_csv("/home/qwang/ss_for_ukb_pqtl/n1384_pqtl_random10K_beta_ukbannot.tsv.gz", sep="\t", compression="gzip")

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
def bin_pip(pip):
    if pip<0.001: return 0
    elif pip<0.01: return 0.001
    elif pip<0.1: return 0.01
    elif pip<0.5: return 0.1
    elif pip<0.9: return 0.5
    elif pip<1: return 0.9
    elif pip==1: return 1
    #elif pip<=1: return 0.9
    else: return np.nan
pj = pd.read_csv("~/Downloads/n1384_pqtl_pip0001_ukbbetaannot.tsv.gz", sep="\t", compression="gzip")
pj["p_min_pip"] = np.minimum(pj.pip_fm, pj.pip_sus)
pj["pqtl_pip_bin"] = pj.p_min_pip.apply(lambda x: bin_pip(x))
pj.pqtl_pip_bin.value_counts() #0 -> other QTL PIPが入っているだけ. can be removed
pj = pj[pj.p_min_pip>=0.001]
pj = pj[~pj.rsid.str.startswith("chrX")] #removing chrX as those are not in UKB

rd = pd.read_csv("~/Downloads/n1384_pqtl_random10K_beta_ukbannot.tsv.gz", sep="\t", compression="gzip")
rd = rd[~rd.variant_hg19.str.startswith("X")] #as a result nothing is removed but just in case
pops = ["eas", "eur", "afr"]
pips = pj.pqtl_pip_bin.sort_values().unique()
colors = ["tab:grey", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
xlabels = ["Random", r"0.001$\leq$PIP<0.01", r"0.01$\leq$PIP<0.1", r"0.1$\leq$PIP<0.5", r"0.5$\leq$PIP<0.9", r"0.9$\leq$PIP<1","PIP=1"]
m = min(pj.slope.min(), pj.BETA_eas.min(), pj.BETA_eur.min(), pj.BETA_afr.min())
M = max(pj.slope.max(), pj.BETA_eas.max(), pj.BETA_eur.max(), pj.BETA_afr.max())
fig, ax = plt.subplots(3, 7, sharex=True, sharey=True, figsize=(12, 5.5))
for i in range(len(pops)):
    for j in range(len(pips)+1):#including randoms
        ax[i,j].axvline(x=0, linestyle="--", color="black", linewidth=0.5, zorder=-2)
        ax[i,j].axhline(y=0, linestyle="--", color="black", linewidth=0.5, zorder=-2)
        ax[i,j].plot([m,M], [m,M], linestyle="--", color="black", linewidth=0.5, zorder=-2)
        ax[i,j].set_aspect("equal")
        if j==0:#random variants
            df = rd
            ax[i,j].set_ylabel("Effect size\nin UKB {0}".format(pops[i].upper()))
        else:
            df = pj[pj.pqtl_pip_bin == pips[j - 1]]
        x, y, yerr, xerr = df.slope, df["BETA_{0}".format(pops[i])], df["SE_{0}".format(pops[i])], df.slope_se
        nonna = ~y.isna()
        x, y, yerr, xerr = x[nonna], y[nonna], yerr[nonna], xerr[nonna]
        ax[i,j].errorbar(x, y, yerr, xerr, fmt="o", color=colors[j], alpha=0.75, mec="black")
        ax[i,j].text(x=m, y=M*0.8, s="r={0}".format(np.round(stats.pearsonr(x,y)[0],3)))
        if i==2:
            ax[i,j].set_xlabel(xlabels[j])
ax[0,0].set_xlim([m*1.05, M*1.05])
ax[0,0].set_ylim([m*1.05, M*1.05])
plt.savefig("/Users/qingbowang/Downloads/vsukb_beta.png", dpi=500)
plt.savefig("/Users/qingbowang/Downloads/vsukb_beta.pdf", dpi=500)
plt.clf()

#Also p-value distribution too
#First need to add "missing" info:
rdexists = pd.read_csv("~/Downloads/n1384_pqtl_random10K_beta_ukbexists.tsv.gz", sep="\t")
rdexists.set_index(["variant_id", "gene_id"], inplace=True)
rd.set_index(["variant_id", "gene_id"], inplace=True)
rd = rd.join(rdexists.iloc[:,-3:])
rd.loc[~rd.is_in_ukb_afr,"LOG10P_afr"] = -1 #as a flag of non-existence
rd.loc[~rd.is_in_ukb_eur,"LOG10P_eur"] = -1 #as a flag of non-existence
rd.loc[~rd.is_in_ukb_eas,"LOG10P_eas"] = -1 #as a flag of non-existence
rd["LOG10P_afr"] = rd["LOG10P_afr"].fillna(0) #as a flag of >0.05
rd["LOG10P_eur"] = rd["LOG10P_eur"].fillna(0) #as a flag of >0.05
rd["LOG10P_eas"] = rd["LOG10P_eas"].fillna(0) #as a flag of >0.05
def bin_log10pval(x):
    if x==-1:
        return -1 #as a flag of missingness
    elif x < -np.log10(0.05):
        return 0 #as a flag of <0.05
    elif x < -np.log10(5e-8):
        return 1 #as a flag of <5e-8
    elif x >= -np.log10(5e-8):
        return 2
    else:
        return (np.nan) #some error flag
rd["LOG10P_afr_bin"] = rd["LOG10P_afr"].apply(lambda x: bin_log10pval(x))
rd["LOG10P_eur_bin"] = rd["LOG10P_eur"].apply(lambda x: bin_log10pval(x))
rd["LOG10P_eas_bin"] = rd["LOG10P_eas"].apply(lambda x: bin_log10pval(x))
rd = rd[~rd.variant_id_hg38.str.startswith("chrX")] #remove X chr (doesn't really exist that much but in case..)
st_rd_afr = rd.LOG10P_afr_bin.value_counts()
st_rd_eur = rd.LOG10P_eur_bin.value_counts()
st_rd_eas = rd.LOG10P_eas_bin.value_counts()

#and the cases:
pjexists = pd.read_csv("~/Downloads/n1384_pqtl_pip0001_ukbexists.tsv.gz", sep="\t")
pjexists.set_index(["rsid"], inplace=True)
pj.set_index(["rsid"], inplace=True)
pj = pj.join(pjexists.iloc[:,-3:])
#remove X chr
pj = pj[~pj.variant_id_hg38.str.startswith("chrX")]
pj.loc[~pj.is_in_ukb_afr.fillna(False),"LOG10P_afr"] = -1 #as a flag of non-existence
pj.loc[~pj.is_in_ukb_eur.fillna(False),"LOG10P_eur"] = -1 #as a flag of non-existence
pj.loc[~pj.is_in_ukb_eas.fillna(False),"LOG10P_eas"] = -1 #as a flag of non-existence
pj["LOG10P_afr"] = pj["LOG10P_afr"].fillna(0) #as a flag of >0.05
pj["LOG10P_eur"] = pj["LOG10P_eur"].fillna(0) #as a flag of >0.05
pj["LOG10P_eas"] = pj["LOG10P_eas"].fillna(0) #as a flag of >0.05
pj["LOG10P_afr_bin"] = pj["LOG10P_afr"].apply(lambda x: bin_log10pval(x))
pj["LOG10P_eur_bin"] = pj["LOG10P_eur"].apply(lambda x: bin_log10pval(x))
pj["LOG10P_eas_bin"] = pj["LOG10P_eas"].apply(lambda x: bin_log10pval(x))

st_case_afr = pj.groupby(["pqtl_pip_bin","LOG10P_afr_bin"]).size().unstack().fillna(0).astype(int)
st_case_eur = pj.groupby(["pqtl_pip_bin","LOG10P_eur_bin"]).size().unstack().fillna(0).astype(int)
st_case_eas = pj.groupby(["pqtl_pip_bin","LOG10P_eas_bin"]).size().unstack().fillna(0).astype(int)

st_rd_afr.name = 0
st_rd_eur.name = 0
st_rd_eas.name = 0
st_afr = pd.concat([st_rd_afr,st_case_afr.T], axis=1).T
st_eur = pd.concat([st_rd_eur,st_case_eur.T], axis=1).T
st_eas = pd.concat([st_rd_eas,st_case_eas.T], axis=1).T

st_afr.to_csv("~/Downloads/ukb_pval_stats_afr.tsv", sep="\t")
st_eur.to_csv("~/Downloads/ukb_pval_stats_eur.tsv", sep="\t")
st_eas.to_csv("~/Downloads/ukb_pval_stats_eas.tsv", sep="\t")

#and plot:
st_afr = pd.read_csv("~/Downloads/ukb_pval_stats_afr.tsv", sep="\t", index_col=0)
st_eur = pd.read_csv("~/Downloads/ukb_pval_stats_eur.tsv", sep="\t", index_col=0)
st_eas = pd.read_csv("~/Downloads/ukb_pval_stats_eas.tsv", sep="\t", index_col=0)

sts = {}
sts["afr"] = st_afr
sts["eur"] = st_eur
sts["eas"] = st_eas
ys = {}
yerrs = {}
for pop in ["afr", "eur", "eas"]:
    ys[pop] = (sts[pop].T / sts[pop].sum(axis=1)).T
    yerrs[pop] = np.sqrt((ys[pop]*(1-ys[pop])).T/ sts[pop].sum(axis=1)).T
    ys[pop].columns.name = "p-value"
    ys[pop].columns = ["Missing", "5e-2<p", "5e-8<p<5e-2", "p<5e-8"]
    #ys[pop] = ys[pop].iloc[:,::-1]

#Plot:
import matplotlib.cm as cm
vir = cm.viridis
cs = ["tab:grey", vir(0.1), vir(0.4), vir(0.75)]
fig, ax = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(10, 6))
ys["eas"].plot.bar(stacked=True, color=cs, ax=ax[0])
ys["eur"].plot.bar(stacked=True, color=cs, ax=ax[1])
ys["afr"].plot.bar(stacked=True, color=cs, ax=ax[2])
ax[0].set_ylabel("Fraction\n(UKB EAS)")
ax[1].set_ylabel("Fraction\n(UKB EUR)")
ax[2].set_ylabel("Fraction\n(UKB AFR)")
ax[2].set_xticks(np.arange(ys["eas"].shape[0]))
ax[2].set_xticklabels(xlabels, rotation=30)
ax[2].set_xlabel("pQTL PIP in our dataset (JCTF)")
ax[1].legend(title="pQTL p-value", loc="lower left", bbox_to_anchor=(1.02,0))
ax[0].legend().set_visible(False)
ax[2].legend().set_visible(False)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/ukb_pval_nain.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/ukb_pval_nain.pdf', bbox_inches='tight', dpi=500)
plt.clf()