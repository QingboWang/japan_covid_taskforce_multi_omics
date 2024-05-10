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
def bin_pip(pip):
    if pip<0.001: return 0
    elif pip<0.01: return 0.001
    elif pip<0.1: return 0.01
    elif pip<0.5: return 0.1
    elif pip<0.9: return 0.5
    elif pip<1: return 0.9
    elif pip==1: return 1
    else: return np.nan

#simply take the enrichment without filtering to any PIP threshold
fm = {}
sus = {}
jn = {}#stands for "joined"
for qtl in ["e","p","es","ps"]:
    fm[qtl] = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n998_{0}qtl_finemap_pip0001.txt".format(qtl), sep=' ')
    sus[qtl] = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n998_{0}qtl_susie_pip0001.txt".format(qtl), sep='\t')
    fm[qtl] = fm[qtl]
    sus[qtl] = sus[qtl]
    fm[qtl]["variant_id_hg38"] = fm[qtl].rsid.str.split("_").str[0]
    fm[qtl]["gene_id"] = fm[qtl].rsid.str.split("_").str[-1]
    fm[qtl].set_index(["variant_id_hg38","gene_id"], inplace=True)
    sus[qtl]["variant_id_hg38"] = sus[qtl].rsid.str.split("_").str[0]
    sus[qtl]["gene_id"] = sus[qtl].rsid.str.split("_").str[-1]
    sus[qtl].set_index(["variant_id_hg38","gene_id"], inplace=True)
    jn[qtl] = fm[qtl].join(sus[qtl].pip, how="inner") #since we will be interested in min pip>0.01 anyways
    jn[qtl].fillna(0, inplace=True) #anything lower than 000001 is 0
    jn[qtl].columns = ["rsid", "{0}_pip_fm".format(qtl), "{0}_pip_sus".format(qtl)]
    jn[qtl]["{0}_min_pip".format(qtl)] = np.minimum(jn[qtl]["{0}_pip_fm".format(qtl)], jn[qtl]["{0}_pip_sus".format(qtl)])
    jn[qtl] = jn[qtl][jn[qtl]["{0}_min_pip".format(qtl)]>=0.001]
    print ("done {0}, {1}".format(qtl, tm.ctime()))


#Unique var-genes:
for qtl in ["e","p","es","ps"]:
    jn[qtl]["Gene"] = jn[qtl].index.get_level_values(1).str.split("\\.").str[0]
    jn[qtl].reset_index(inplace=True)
    jn[qtl].set_index(["variant_id_hg38", "Gene"], inplace=True)
univg = jn["e"].index.union(jn["p"].index.union(jn["es"].index.union(jn["ps"].index)))
#vep-annotation:
consq = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/vepannot_multiqtl_pip0001_consequence_pervg.tsv",sep='\t', compression="gzip", index_col=[0,1])
for qtl in ["e","p","es","ps"]:
    jn[qtl] = jn[qtl].join(consq.Consequence, how="left")
#Also annotate the list of TF binding variants, per variant annotation:
univ = univg.get_level_values(0)
vps = []
cols_to_read = ["#Uploaded_variation" , "Consequence"]
for chr in list(range(1,23))+["X"]:
    vp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/all_n1419_vepped/n1419_vepped_chr{0}.txt".format(chr), sep='\t', usecols=cols_to_read)
    vp = vp[vp.Consequence.str.contains("TF_binding")]
    vp["variant_id_hg38"] = vp["#Uploaded_variation"].str.split("_").str[0]
    vp.set_index(["variant_id_hg38"], inplace=True)
    intsct = vp.index.intersection(univ)
    vp = vp.loc[intsct,:]
    vps = vps + list(vp.index.unique())
    if chr==1:
        print (vp)
    print ("done chr{0}, {1}".format(chr, tm.ctime()))
vps = pd.Series(vps)
vps.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/multiqtl_pip0001_tfvars.tsv.gz", sep='\t', compression='gzip', index=False)
#annotate:
vps = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/multiqtl_pip0001_tfvars.tsv.gz", sep='\t',
                  compression='gzip', index_col=0)
vps["TF_binding"] = True
vps.index.names = ["variant_id_hg38"]
for qtl in ["e","p","es","ps"]:
    jn[qtl] = jn[qtl].join(vps.TF_binding, how="left", on="variant_id_hg38")

#Filtering out CLPP>0.01 to make thing clearer
ep = jn["e"].join(jn["p"].p_min_pip, how="inner")
(ep.e_min_pip*ep.p_min_pip).sort_values().tail(40)
tofilt = ep[(ep.e_min_pip*ep.p_min_pip)>0.01].index #0.01 it is!
#and get the enrichment:
annots = ["intron_variant", "upstream_gene_variant", "downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"]
annots += ["synonymous_variant", "missense_variant", "inframe_deletion|inframe_insertion|start_lost|stop_lost", "frameshift_variant|stop_gained"]
ns = {}
Ns = {}
frac = {}
err = {}
for qtl in ["e","p","es","ps"]:
    ns[qtl] = {}
    Ns[qtl] = {}
    frac[qtl] = {}
    err[qtl] = {}
    jn[qtl]["pip_bin"] = jn[qtl]["{0}_min_pip".format(qtl)].apply(lambda x: bin_pip(x))
    for annot in annots:
        jn[qtl][annot] = jn[qtl].Consequence.fillna("NA").str.contains(annot)
        #tb = jn[qtl].groupby(["pip_bin", annot]).size().unstack().fillna(0).astype(int)
        tb = jn[qtl].loc[jn[qtl].index.difference(tofilt),:].groupby(["pip_bin", annot]).size().unstack().fillna(0).astype(int) #Filtering by CLPP<0.01
        f = tb[True] / tb.sum(axis=1)
        e = np.sqrt(f*(1-f)/tb.sum(axis=1))
        ns[qtl][annot] = tb[True]
        Ns[qtl][annot] = tb.sum(axis=1)
        frac[qtl][annot] = f
        err[qtl][annot] = e
    frac[qtl] = pd.DataFrame(frac[qtl])
    err[qtl] = pd.DataFrame(err[qtl])
#For TF binding as well:
for qtl in ["e","p","es","ps"]:
    jn[qtl]["TF_binding"] = jn[qtl]["TF_binding"].fillna(False)
    tb = jn[qtl].loc[jn[qtl].index.difference(tofilt),:].groupby(["pip_bin", "TF_binding"]).size().unstack().fillna(0).astype(int) #Filtering by CLPP<0.01
    f = tb[True] / tb.sum(axis=1)
    e = np.sqrt(f*(1-f)/tb.sum(axis=1))
    ns[qtl]["TF_binding"] = tb[True]
    Ns[qtl]["TF_binding"] = tb.sum(axis=1)
    frac[qtl]["TF_binding"] = f
    err[qtl]["TF_binding"] = e
#save:
for qtl in ["e","p","es","ps"]:
    ns[qtl] = pd.DataFrame(ns[qtl])
    Ns[qtl] = pd.DataFrame(Ns[qtl])
    frac[qtl] = pd.DataFrame(frac[qtl])
    err[qtl] = pd.DataFrame(err[qtl])
    ns[qtl].to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/nnumer_allqtl_funcannot_{0}.tsv".format(qtl), sep="\t")
    Ns[qtl].to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ndenom_allqtl_funcannot_{0}.tsv".format(qtl), sep="\t")
    frac[qtl].to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/frac_allqtl_funcannot_{0}.tsv".format(qtl), sep="\t")
    err[qtl].to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/err_allqtl_funcannot_{0}.tsv".format(qtl), sep="\t")

#Get the numbers in the background (when no stratification by PIP)
for chr in list(range(1,23))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/esqtl_sumstats/esqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
    es = pd.read_csv(input_file, sep="\t")
    es["variant_id_hg38"] = es.variant_id.str.split("_").str[0]
    es["Gene"] = es.gene_id.str.split("\\.").str[0]
    es.set_index(["variant_id_hg38","Gene"], inplace=True)
    vp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/all_n1419_vepped/n1419_vepped_chr{0}.txt".format(chr),sep='\t')
    vp["variant_id_hg38"] = vp["#Uploaded_variation"].str.split("_").str[0]
    vp.set_index(["variant_id_hg38","Gene"], inplace=True)
    cons = vp.groupby(["variant_id_hg38", "Gene"]).Consequence.apply(lambda x: ','.join(x))
    es = es.join(cons, how="left")
    es.Consequence.fillna("NA", inplace=True)
    uniqueannot = pd.Series((",".join(es.Consequence.unique())).split(",")).unique()
    tb = {}
    for annot in uniqueannot:
        t = es.Consequence.str.contains(annot).value_counts()
        tb[annot] = t
    es["inframe"] = es.Consequence.str.contains("inframe_deletion|inframe_insertion|start_lost|stop_lost")
    es["ptv"] = es.Consequence.str.contains("frameshift_variant|stop_gained")
    for annot in ["inframe","ptv"]:
        t = es[annot].value_counts()
        tb[annot] = t
    tb = pd.DataFrame(tb)
    tb.fillna(0).astype(int).to_csv("/Users/qingbowang/Desktop/tmp/vepannot_esqtl_baseline_stats_chr{0}_updated.tsv".format(chr), sep='\t')
    print ("done chr{0}, {1}".format(chr, tm.ctime()))

tb = pd.read_csv("/Users/qingbowang/Desktop/tmp/vepannot_esqtl_baseline_stats_chr{0}_updated.tsv".format(1), sep='\t', index_col=0)
for chr in list(range(2,23))+["X"]:
    tb = tb.add(pd.read_csv("/Users/qingbowang/Desktop/tmp/vepannot_esqtl_baseline_stats_chr{0}_updated.tsv".format(chr), sep='\t', index_col=0), fill_value=0)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_esqtl_baseline_stats_updated.tsv".format(chr), sep='\t')
#numerator:
vps = []
cols_to_read = ["#Uploaded_variation" , "Consequence"]
for chr in list(range(1,23))+["X"]:
    vp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/all_n1419_vepped/n1419_vepped_chr{0}.txt".format(chr), sep='\t', usecols=cols_to_read)
    vp = vp[vp.Consequence.str.contains("TF_binding")]
    vp["variant_id_hg38"] = vp["#Uploaded_variation"].str.split("_").str[0]
    vps = vps + list(vp.variant_id_hg38.unique())
    print ("done chr{0}, {1}".format(chr, tm.ctime()))
vps = pd.Series(vps)
vps.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/multiqtl_all_tfvars.tsv.gz", sep='\t', compression='gzip', index=False)
vps = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/multiqtl_all_tfvars.tsv.gz",
                  sep='\t', compression='gzip', index_col=0)
vps["TF_binding"] = True
vps.index.names = ["variant_id_hg38"]
N0 = 0
n = 0
for chr in list(range(1,23))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/esqtl_sumstats/esqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr)
    df = pd.read_csv(input_file, sep="\t", usecols=["variant_id"])
    df.index = df.variant_id.str.split("_").str[0]
    df.index.names = ["variant_id_hg38"]
    N0 += df.shape[0]
    df = df.join(vps, how="inner")
    n += df.shape[0]
    print ("done chr {0}, so far N0={1}, n={2}, {3}".format(chr, N0, n, tm.ctime()))
bg = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_esqtl_baseline_stats_updated.tsv", sep='\t', index_col=0)
bg["TF_binding"] = [N0-n, n]
bg.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_esqtl_baseline_stats_updated_wtf.tsv", sep='\t')

#Plot enrichment:
bg = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_esqtl_baseline_stats_updated_wtf.tsv", sep='\t', index_col=0)
qtl = "es"
ns = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/nnumer_allqtl_funcannot_{0}.tsv".format(qtl),sep="\t", index_col=0)
Ns = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ndenom_allqtl_funcannot_{0}.tsv".format(qtl),sep="\t", index_col=0)

ns.rename(columns = {'inframe_deletion|inframe_insertion|start_lost|stop_lost':"inframe", 'frameshift_variant|stop_gained':"ptv"}, inplace=True)
Ns.rename(columns = {'inframe_deletion|inframe_insertion|start_lost|stop_lost':"inframe", 'frameshift_variant|stop_gained':"ptv"}, inplace=True)
ns.loc["all"] = bg.loc[True,ns.columns]
Ns.loc["all"] = bg[Ns.columns].sum(axis=0)
ns.loc[0] = ns.loc["all",:]-ns.iloc[:-1,:].sum()
Ns.loc[0] = Ns.loc["all",:]-Ns.iloc[:-1,:].sum()
ns = ns.iloc[[6,7,0,1,2,3,4,5],:]
Ns = Ns.iloc[[6,7,0,1,2,3,4,5],:]
p_pip_given_annot = ns/ns.iloc[0,:] #This is p(PIP | annot)
p_pip_random = Ns.iloc[:,0]/Ns.iloc[0,0]
enr_es_correct = (p_pip_given_annot.T/p_pip_random).T #This is p(PIP | annot) / p(PIP)
frac = ns/Ns #This is the p(annot | PIP)
err = np.sqrt((frac*(1-frac)/Ns))
enr_es = frac.iloc[1:,:]/frac.iloc[0,:] #This is p(annot | PIP) / p(annot)
enr_err_es = err.iloc[1:,:]/frac.iloc[0,:]
qtl = "ps"
ns = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/nnumer_allqtl_funcannot_{0}.tsv".format(qtl),sep="\t", index_col=0)
Ns = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ndenom_allqtl_funcannot_{0}.tsv".format(qtl),sep="\t", index_col=0)
ns.rename(columns = {'inframe_deletion|inframe_insertion|start_lost|stop_lost':"inframe", 'frameshift_variant|stop_gained':"ptv"}, inplace=True)
Ns.rename(columns = {'inframe_deletion|inframe_insertion|start_lost|stop_lost':"inframe", 'frameshift_variant|stop_gained':"ptv"}, inplace=True)
ns.loc["all"] = bg.loc[True,ns.columns]
Ns.loc["all"] = bg[Ns.columns].sum(axis=0)
ns.loc[0] = ns.loc["all",:]-ns.iloc[:-1,:].sum()
Ns.loc[0] = Ns.loc["all",:]-Ns.iloc[:-1,:].sum()
ns = ns.iloc[[6,7,0,1,2,3,4,5],:]
Ns = Ns.iloc[[6,7,0,1,2,3,4,5],:]
frac = ns/Ns
err = np.sqrt((frac*(1-frac)/Ns))
enr_ps = frac.iloc[1:,:]/frac.iloc[0,:]
enr_err_ps = err.iloc[1:,:]/frac.iloc[0,:]

#Plot:
cols_use = ["TF_binding", "5_prime_UTR_variant", "3_prime_UTR_variant", 'ptv', 'missense_variant']
idx = ["TF binding", "5'UTR variant", "3'UTR variant", 'protein-truncating\nvariant','missense\nvariant']
c7 = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tk = ["[0,0.001)","[0.001,0.01)", "[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1]", "1"]
ys_use_es = enr_es.loc[:,cols_use].T
yerrs_use_es = enr_err_es.loc[:,cols_use].T
ys_use_ps = enr_ps.loc[:,cols_use].T
yerrs_use_ps = enr_err_ps.loc[:,cols_use].T
x = np.arange(ys_use_es.shape[0])
fig = plt.figure(figsize=(8,2.7))
i = 0
for pip in ys_use_es.columns:
    y = ys_use_es[pip]
    yerr = yerrs_use_es[pip]
    plt.errorbar(x-0.45+0.15*i-0.03, y, yerr, fmt="^", color=c7[i], mec='black')
    y = ys_use_ps[pip]
    yerr = yerrs_use_ps[pip]
    plt.errorbar(x - 0.45 + 0.15 * i + 0.03, y, yerr, fmt="v", color=c7[i], mec='black')
    i += 1
#overlay lines:
j = 0
for categ in ys_use_es.index:
    plt.plot(j-0.45 - 0.03 + np.arange(ys_use_es.shape[1])/10*1.5, ys_use_es.loc[categ,:], linestyle="--", linewidth=0.5, color="tab:cyan", zorder=-2)
    plt.plot(j-0.45 + 0.03 + np.arange(ys_use_ps.shape[1]) / 10*1.5, ys_use_ps.loc[categ, :], linestyle="--", linewidth=0.5, color="tab:olive", zorder=-2)
    j += 1
#just for legend:
i = 0
for pip in ys_use_es.columns:
    plt.errorbar(0, 100**-10, 0, fmt="o", color=c7[i], mec='black', label=tk[i])
    i += 1
color_legend = plt.legend(title='QTL PIP bin',loc='upper left', bbox_to_anchor=(1.01,1.1))
plt.gca().add_artist(color_legend)
# Create a list of handles for the lines that are plotted with plt.plot()
shape_handles = []
shape_handles.append(plt.plot(0,10**-10, marker="^", linestyle="--", linewidth=0.5, color="tab:cyan", mec='black', markerfacecolor='black', label="mRNA specific")[0])
shape_handles.append(plt.plot(0,100**-10, marker="v", linestyle="--", linewidth=0.5, color="tab:olive", mec='black', markerfacecolor='black', label="Protein specific")[0])
# Pass the list of handles to the shape_legend
shape_legend = plt.legend(handles=shape_handles,title='QTL classification',loc='upper left', bbox_to_anchor=(1.01,0.1 ))
plt.gca().add_artist(shape_legend)
plt.xticks(x, idx, rotation=30)
plt.xlabel("Annotation")
plt.ylabel("Enrichment")
plt.yscale('log')
plt.ylim([10**-0.75 ,10**4.1 ])
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/espsqtl_vep_enr_updated.png', bbox_inches='tight', dpi=500)
plt.clf()

#CLPP: sort by different PIPs and output as tables:
cols_use = ["TF_binding", "5_prime_UTR_variant", "3_prime_UTR_variant", 'ptv', 'missense_variant']
for qtl in ["e","p","es","ps"]:
    jn[qtl].rename(columns={'inframe_deletion|inframe_insertion|start_lost|stop_lost': "inframe",
                       'frameshift_variant|stop_gained': "ptv"}, inplace=True)
    jn[qtl].sort_values(by="{0}_min_pip".format(qtl), ascending=False, inplace=True)
    bests = jn[qtl].iloc[:2001, :]
    bests["N"] = np.arange(bests.shape[0])+1
    for col in cols_use:
        bests["{0}_cumsum".format(col)] = bests[col].cumsum()
        bests["{0}_frac".format(col)] = bests["{0}_cumsum".format(col)]/bests.N
        bests["{0}_err".format(col)] = np.sqrt(bests["{0}_frac".format(col)]*(1-bests["{0}_frac".format(col)])/bests.N)
    bests.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/{0}qtl_funcannot_cumfrac.tsv.gz".format(qtl), sep='\t', compression='gzip')
ep = jn["e"].join(jn["p"].p_min_pip, how="inner")
ep["clpp"] = ep.e_min_pip*ep.p_min_pip
ep.rename(columns={'inframe_deletion|inframe_insertion|start_lost|stop_lost': "inframe",
                        'frameshift_variant|stop_gained': "ptv"}, inplace=True)
ep.sort_values(by="clpp", ascending=False, inplace=True)
bests = ep.iloc[:2001, :]
bests["N"] = np.arange(bests.shape[0]) + 1
for col in cols_use:
    bests["{0}_cumsum".format(col)] = bests[col].cumsum()
    bests["{0}_frac".format(col)] = bests["{0}_cumsum".format(col)] / bests.N
    bests["{0}_err".format(col)] = np.sqrt(
        bests["{0}_frac".format(col)] * (1 - bests["{0}_frac".format(col)]) / bests.N)
bests.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/clpp_funcannot_cumfrac.tsv.gz".format(qtl), sep='\t', compression='gzip')
#agg them all and plot (for this, we don't even need backgrounds...):
bests = {}
for qtl in ["eqtl","pqtl","esqtl","psqtl", "clpp"]:
    bests[qtl] = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/{0}_funcannot_cumfrac.tsv.gz".format(qtl),sep='\t', nrows=10000) #10000 is enough
funcnames = ["TF binding", "5'UTR variant", "3'UTR variant", 'protein-truncating variant','missense variant']
qtlnames = ["mRNA", "Protein", "mRNA specific", "Protein specific", "Colocalizing"]
colors = ["tab:blue","tab:green","tab:cyan", "tab:olive", "tab:orange"]
fig, ax = plt.subplots(3, 2, figsize = (8,6))
cnt = 0 #total cnt
r = 0 #row
c = 0 #column
for col in cols_use:
    i = 0
    for qtl in ["eqtl","pqtl","esqtl","psqtl", "clpp"]:
        x = bests[qtl].N[50:1929]
        y = bests[qtl]["{0}_frac".format(col)][50:1929]
        yerr = bests[qtl]["{0}_err".format(col)][50:1929]
        if (r == 2) and (c == 0):
            ax[r,c].plot(x, y, label=qtlnames[i], color=colors[i])
        else:
            ax[r, c].plot(x, y, color=colors[i]) #for plotting
        ax[r,c].fill_between(x, y - yerr, y + yerr, alpha=0.2, color=colors[i])
        i += 1
    ax[r, c].text(1000, ax[r,c].get_ylim()[1]*0.95, funcnames[cnt], horizontalalignment='center',verticalalignment='center')
    ax[r, c].spines.right.set_visible(False)
    ax[r, c].spines.top.set_visible(False)
    ax[r, c].set_ylabel("Fraction")
    ax[r, c].invert_xaxis()
    if r<2: r += 1
    else:
        r = 0
        c += 1
    cnt += 1
ax[1,0].set_xlabel("Best x variant-genes")
ax[1,1].set_xlabel("Best x variant-genes")
fig.delaxes(ax[2, 1])
plt.tight_layout()
plt.legend(title="QTL Classification", bbox_to_anchor=(1.2, 1.0))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/best_qtls_funcenr.png', bbox_inches='tight', dpi=500)
plt.clf()


