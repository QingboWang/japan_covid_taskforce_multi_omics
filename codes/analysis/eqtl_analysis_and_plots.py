
import pandas as pd
import numpy as np
import time as tm

##Comparison with previous release:

#n465 (previous release):
minp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_rna_releasedata_cis_eqtls.tsv.gz", sep='\t', compression='gzip')
minp_n465 = minp.groupby("gene_id").pval_nominal.min() #only if it were p<0.05, which is fine in practice
minp_n465.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n465.tsv", sep="\t")
#n1019, old:
minps_old = []
for chr in range(1,23):
    df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/eqtl_sumstats/n1019.sevon.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t')
    minps_old.append(df.groupby("gene_id").pval_nominal.min())
    print ("done {0}, {1}".format(chr, tm.ctime()))
minps_old = pd.concat(minps_old)
minps_old.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019.tsv", sep="\t")
#n1019, new
minps = []
for chr in list(range(1,23))+["X"]:
    df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t')
    minps.append(df.groupby("gene_id").pval_nominal.min())
    print ("done {0}, {1}".format(chr, tm.ctime()))
minps = pd.concat(minps)
minps.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019_updated.tsv", sep="\t")

#read back and join:
minp_n465 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n465.tsv", sep="\t", index_col=0)
minps_old = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019.tsv", sep="\t", index_col=0)
minps = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019_updated.tsv", sep="\t", index_col=0)
minp_n465.columns = ["n465"]
minps_old.columns = ["n1019_old"]
minps.columns = ["n1019"]
minps = minps.join(pd.DataFrame(minp_n465), how="left").join(pd.DataFrame(minps_old), how="left")
#save:
minps.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019_vs_n465.tsv", sep="\t")

#comparison table:
def bin_pval(p):
    if p<10**-100: return (100)
    elif p<10**-50: return (50)
    elif p<10**-30: return (30)
    elif p<10**-20: return (20)
    elif p<10**-10: return (10)
    elif p<10**-9: return (9)
    elif p<10**-8: return (8)
    elif p<10**-7: return (7)
    elif p<10**-6: return (6)
    elif p<10**-5: return (5)
    elif p<10**-4: return (4)
    elif p<10**-3: return (3)
    elif p<10**-2: return (2)
    else: return (0)
minps = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019_vs_n465.tsv", sep="\t", index_col=0)
minps = minps[minps.isna().sum(axis=1)==0] #filtering NAs
minps["n1019_bin"] = minps.n1019.apply(lambda x: bin_pval(x))
minps["n1019_old_bin"] = minps.n1019_old.apply(lambda x: bin_pval(x))
minps["n465_bin"] = minps.n465.apply(lambda x: bin_pval(x))
tb = minps.groupby(["n1019_bin", "n465_bin"]).size().unstack().fillna(0).astype(int)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019_vs_n465_stats.tsv", sep="\t")
tb = minps.groupby(["n1019_bin", "n1019_old_bin"]).size().unstack().fillna(0).astype(int)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019_vs_preimputation_stats.tsv", sep="\t")


##comparison in terms of PIPs

#read the data:
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=" ")
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep="\t")
e_fm["variant_id_hg38"] = e_fm.rsid.str.split("_").str[0]
e_fm["gene_id"] = e_fm.rsid.str.split("_").str[-1]
e_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
e_sus["variant_id_hg38"] = e_sus.rsid.str.split("_").str[0]
e_sus["gene_id"] = e_sus.rsid.str.split("_").str[-1]
e_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)

#check 1: SuSiE vs FM for each.
ej = e_fm.join(e_sus.pip, how="outer")
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej.columns = ["rsid", "pip_fm", "pip_sus"]
def bin_pip(pip):
    if pip<0.001: return 0
    elif pip<0.01: return 0.001
    elif pip<0.1: return 0.01
    elif pip<0.5: return 0.1
    elif pip<0.9: return 0.5
    elif pip<1: return 0.9
    elif pip==1: return 1
    else: return np.nan
ej["fm_bin"] = ej.pip_fm.apply(lambda x: bin_pip(x))
ej["sus_bin"] = ej.pip_sus.apply(lambda x: bin_pip(x))
ej["min_pip"] = np.minimum(ej.pip_fm, ej.pip_sus)
ej["min_pip_bin"] = ej.min_pip.apply(lambda x: bin_pip(x))
#table:
tb = ej.groupby(["fm_bin", "sus_bin"]).size().unstack().fillna(0).astype(int)

#also add number for 0,0:
e_fm_low = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_n_below_pip0001.txt", sep='\s+', index_col=0)
e_sus_low = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_n_below_pip0001.txt", sep='\s+', index_col=0)
tb.iloc[0,0] = e_fm_low.sum() - tb.sum(axis=1)[0] #all the fm bin==0 - (fm bin==0 and sus bin!=0) = (fm_bin==0 and sus_bin==0)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/eqtl_sus_vs_fm_bins.tsv", sep='\t')

#Get the list of variants missing in the previous version
v465 = []
for chk in range(26):
    df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot.txt.gz".format(chk), sep='\t')
    v465 = v465+ list(df.variant_id.unique())
    print ("done chk{0}, {1}".format(chk, tm.ctime()))
pd.Series(v465).to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/unique_variants_n465.tsv", sep="\t", index=False)
v1019 = []
for chr in list(range(1,23))+["X"]:
    df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t')
    v1019 = v1019 + list(df.variant_id.unique())
    print("done chr{0}, {1}".format(chr, tm.ctime()))
pd.Series(v1019).to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/unique_variants_n1019.tsv", sep="\t", index=False)
#Get the list of variants in our analysis that are missing in n465 analysis:
v465 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/unique_variants_n465.tsv", sep="\t")
v465 = pd.DataFrame(v465.unique())
v465.index = v465.iloc[:,0]
v465["exists_in_n465"] = True
v1019 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/unique_variants_n1019.tsv", sep="\t")
v1019 = pd.DataFrame(v1019)
v1019.columns = ["orig"]
v1019.index = v1019.orig.str.replace("_",":")
v1019 = v1019.join(v465.exists_in_n465, how="left")
v1019.index = v1019.orig
v1019.exists_in_n465.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/unique_variants_n1019_exists_in_n465_or_not.tsv", sep="\t")

#Use this to get the stats:
#first get the variant exists or not, p<5e-8 or not flag
vs = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/unique_variants_n1019_exists_in_n465_or_not.tsv", sep="\t", index_col=0)
vs.index.names = ["variant_id"]
minp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_rna_releasedata_cis_eqtls.tsv.gz", sep='\t', compression='gzip')
minp_n465 = minp.groupby("gene_id").pval_nominal.min() #only if it were p<0.05, which is fine in practice
#and also the prev PIP
pips_n465 = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/susie_and_fm_n465_imputed_all_v.tsv.zip", sep='\t')
pips_n465[(pips_n465.pip_fm >= 0.001)&(pips_n465.pip_susie >= 0.001)].to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/susie_and_fm_n465_imputed_min0001.tsv.gz", sep='\t', index=False, compression="gzip") #<0.001 -> don't care
pips_n465 = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/susie_and_fm_n465_imputed_min0001.tsv.gz", sep='\t', index_col=[0,1])
pips_n465.index.names = ["variant_id_dots", "gene_id"]
pips_n465["pip_n465"] = np.minimum(pips_n465.pip_fm, pips_n465.pip_susie)

#The vg list for our current fine-mapping with n1019 re-imputed
vg1019 = []
for chr in list(range(1,23))+["X"]:
    df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t', index_col=[0,1])
    vg1019.append(df)
    print ("done {0}, {1}".format(chr, tm.ctime()))
vg1019 = pd.concat(vg1019)
#annoate eGene info from n465
vg1019 = vg1019.join(minp_n465, how="left", on="gene_id", rsuffix="_min_n465")
print ("done {0}, {1}".format("eGene info n465", tm.ctime()))
#annotate whether variant exists in n465
vg1019 = vg1019.join(vs, how="left", on="variant_id")
print ("done {0}, {1}".format("variant existing info n465", tm.ctime()))
#annotate pips info from n465 (<0.0001 => 0) (str replace takes time..)
vg1019["variant_id_dots"] = vg1019.index.get_level_values(1).str.replace("_",":")
vg1019.set_index("variant_id_dots", append=True, inplace=True)
vg1019.reset_index(level="variant_id",  inplace=True)
vg1019 = vg1019.join(pips_n465.pip_n465, how="left", on=["variant_id_dots", "gene_id"])
print ("done {0}, {1}".format("pip info n465", tm.ctime()))
#annotate pips info from our n1019 (<0.001 => 0)
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=" ")
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep="\t")
e_fm["variant_id_dots"] = e_fm.rsid.str.split("_").str[0]+":"+e_fm.rsid.str.split("_").str[1]
e_fm["gene_id"] = e_fm.rsid.str.split("_").str[-1]
e_fm.set_index(["variant_id_dots","gene_id"], inplace=True)
e_sus["variant_id_dots"] = e_sus.rsid.str.split("_").str[0]+":"+e_sus.rsid.str.split("_").str[1]
e_sus["gene_id"] = e_sus.rsid.str.split("_").str[-1]
e_sus.set_index(["variant_id_dots","gene_id"], inplace=True)
ej = e_fm.join(e_sus.pip, how="inner")
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej.columns = ["rsid", "pip_fm", "pip_sus"]
ej["pip_min"] = np.minimum(ej.pip_fm, ej.pip_sus)
vg1019 = vg1019.join(ej.pip_min, how="left", on=["variant_id_dots", "gene_id"])
print ("done {0}, {1}".format("pip info n1019", tm.ctime()))
vg1019.reset_index(level="variant_id_dots",  inplace=True, drop=True)
vg1019['pval_nominal_min_n465'].fillna(1, inplace=True)
vg1019['pip_n465'].fillna(-1, inplace=True)
vg1019['pip_fm'].fillna(0, inplace=True)
vg1019['pip_sus'].fillna(0, inplace=True)
vg1019['pip_min'].fillna(0, inplace=True)
vg1019['exists_in_n465'].fillna(False, inplace=True)
#and then get the stats:
#for each n1019 PIP bin, classify as [was missing], [exists but was not eGene], [0,0.0001), ...
vg1019["n465_class"] = 100 #label for "haven't annotated"
vg1019.loc[~vg1019.exists_in_n465, "n465_class"] = -2 #label for "was missing"
vg1019.loc[vg1019.exists_in_n465 & (vg1019.pval_nominal_min_n465>5*10**-8), "n465_class"] = -1 #label for "was not eGene"
vg1019.loc[vg1019.exists_in_n465 & (vg1019.pip_n465<0.001) & (vg1019.pval_nominal_min_n465 <= 5*10**-8), "n465_class"] = 0 #labeling the lower
vg1019.loc[vg1019.exists_in_n465 & (vg1019.pip_n465>0.001), "n465_class"] = 0.001 #labeling the lower
vg1019.loc[vg1019.exists_in_n465 & (vg1019.pip_n465>0.01), "n465_class"] = .01 #labeling the lower
vg1019.loc[vg1019.exists_in_n465 & (vg1019.pip_n465>0.1), "n465_class"] = .1 #labeling the lower
vg1019.loc[vg1019.exists_in_n465 & (vg1019.pip_n465>0.5), "n465_class"] = .5 #labeling the lower
vg1019.loc[vg1019.exists_in_n465 & (vg1019.pip_n465>0.9), "n465_class"] = .9 #labeling the lower
vg1019.loc[vg1019.exists_in_n465 & (vg1019.pip_n465==1), "n465_class"] = 1 #labeling the lower
vg1019.n465_class.value_counts()
#Save:
vg1019.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_allvariants_vs_n465.tsv.gz", sep="\t")
vg1019 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_allvariants_vs_n465.tsv.gz", sep="\t")#takes time to read
#and get the sum table:
tb = []
vg1019["pip_n1019"] = vg1019.pip_min #just for name consistency
dfsub = vg1019
for i in np.arange(101) / 100:
    dfsub = dfsub[dfsub.pip_n1019 >= i]
    tb.append(dfsub.n465_class.value_counts())
    print ("done {0}, {1}".format(i, tm.ctime()))
tb = pd.DataFrame(tb)
tb.index = range(101)
tb.columns = ["Variant missing", "Gene not an eGene", "(0,0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tb = tb.fillna(0).astype(int)
print(tb)
print(tb.apply(lambda x: x / sum(x), axis=1))
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/vs_n465_pip_updated.tsv", sep='\t') #updatedはやはりつけた方が良い

#and simple count (bin) for both
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=" ", index_col=0)
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep="\t", index_col=0)
ej = e_fm.join(e_sus.pip, how="inner")
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej.columns = ["pip_fm", "pip_sus"]
ej["pip_min"] = np.minimum(ej.pip_fm, ej.pip_sus)

ejold = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/susie_and_fm_n465_imputed_all_v.tsv.zip", sep='\t')#takes some time to read
ejold["pip_min"] = np.minimum(ejold.pip_fm, ejold.pip_susie)
ej["min_pip_bin"] = ej.pip_min.apply(lambda x: bin_pip(x))
ejold["min_pip_bin"] = ejold.pip_min.apply(lambda x: bin_pip(x))
ejold = ejold[ejold.min_pip_bin>=0.001] #filtering out meaningless ones
ej.min_pip_bin.value_counts().sort_index()
ejold.min_pip_bin.value_counts().sort_index()
st = pd.concat([ej.min_pip_bin.value_counts().sort_index(), ejold.min_pip_bin.value_counts().sort_index()], axis=1)
st.columns = ["n1019", "n465"]
st.index = ["[0.001,0.01)","[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
st.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/vs_n465_pip_number.tsv", sep='\t')

##Also functional validation that current version is as good, by TSS distance
abs_log_tss_distance = vg1019.tss_distance.replace(0,1).apply(lambda x: np.log10(abs(x)))
vg1019["tss_distance_bin"] = 100 #to begin with, for error catch
vg1019.loc[abs_log_tss_distance <= 6, "tss_distance_bin"] = 6 #upper bound
vg1019.loc[abs_log_tss_distance <= 5.5, "tss_distance_bin"] = 5.5
vg1019.loc[abs_log_tss_distance <= 5, "tss_distance_bin"] = 5
vg1019.loc[abs_log_tss_distance <= 4.5, "tss_distance_bin"] = 4.5
vg1019.loc[abs_log_tss_distance <= 4, "tss_distance_bin"] = 4
vg1019.loc[abs_log_tss_distance <= 3.5, "tss_distance_bin"] = 3.5
vg1019.loc[abs_log_tss_distance <= 3, "tss_distance_bin"] = 3
vg1019.loc[abs_log_tss_distance <= 2.5, "tss_distance_bin"] = 2.5
vg1019.loc[abs_log_tss_distance <= 2, "tss_distance_bin"] = 2

vg1019["pip_n1019_bin"] = 0
vg1019.loc[vg1019.pip_n1019>0.001, "pip_n1019_bin"] = .001 #labeling the lower
vg1019.loc[vg1019.pip_n1019>0.01, "pip_n1019_bin"] = .01 #labeling the lower
vg1019.loc[vg1019.pip_n1019>0.1, "pip_n1019_bin"] = .1 #labeling the lower
vg1019.loc[vg1019.pip_n1019>0.5, "pip_n1019_bin"] = .5 #labeling the lower
vg1019.loc[vg1019.pip_n1019>0.9, "pip_n1019_bin"] = .9 #labeling the lower
vg1019.loc[vg1019.pip_n1019==1, "pip_n1019_bin"] = 1 #labeling the lower

tb = vg1019.groupby(["tss_distance_bin", "pip_n1019_bin"]).size().unstack().fillna(0).astype(int)
tb_called = vg1019[vg1019.exists_in_n465&(vg1019.pval_nominal_min_n465<5*10**-8)].groupby(["tss_distance_bin", "pip_n1019_bin"]).size().unstack().fillna(0).astype(int)
tb_was_missing = vg1019[~vg1019.exists_in_n465].groupby(["tss_distance_bin", "pip_n1019_bin"]).size().unstack().fillna(0).astype(int)
tb_wasnt_egene = vg1019[(vg1019.pval_nominal_min_n465>5*10**-8)].groupby(["tss_distance_bin", "pip_n1019_bin"]).size().unstack().fillna(0).astype(int)
tb_was_missing_egene = vg1019[(~vg1019.exists_in_n465)&(vg1019.pval_nominal_min_n465<5*10**-8)].groupby(["tss_distance_bin", "pip_n1019_bin"]).size().unstack().fillna(0).astype(int)
tb_wasnt_egene_exists = vg1019[vg1019.exists_in_n465 & (vg1019.pval_nominal_min_n465>5*10**-8)].groupby(["tss_distance_bin", "pip_n1019_bin"]).size().unstack().fillna(0).astype(int) #Existsに絞っている.

#save these
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_updated.tsv", sep='\t')
tb_called.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_intersection_updated.tsv", sep='\t')
tb_was_missing.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_prev_missingvars_updated.tsv", sep='\t')
tb_wasnt_egene.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_prev_nonegenes_updated.tsv", sep='\t')
tb_was_missing_egene.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_prev_missingvars_egene_updated.tsv", sep='\t')
tb_wasnt_egene_exists.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_prev_nonegenes_existingvars_updated.tsv", sep='\t')

#Also go EMS distribution:
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=" ", index_col=0)
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep="\t", index_col=0)
ej = e_fm.join(e_sus.pip, how="inner")
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej.columns = ["pip_fm", "pip_sus"]
ej["pip_min"] = np.minimum(ej.pip_fm, ej.pip_sus)
ejold = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/susie_and_fm_n465_imputed_all_v.tsv.zip", sep='\t')#takes some time to read
ejold["pip_min"] = np.minimum(ejold.pip_fm, ejold.pip_susie)
ej["min_pip_bin"] = ej.pip_min.apply(lambda x: bin_pip(x))
ejold["min_pip_bin"] = ejold.pip_min.apply(lambda x: bin_pip(x))
ejold = ejold[ejold.min_pip_bin>=0.001] #filtering out meaningless ones
#GTEx PIP and EMS
gpip = pd.read_csv("/Users/qingbowang/Desktop/resources/Whole_Blood_allpairs_gtex_pip.tsv.gz", sep="\t") #takes some time, large file
gpip["gene"] = gpip.gene_id.str.split("\\.").str[0] #takes some time..
gpip.set_index(["variant_id", "gene"], inplace=True)
ej["variant_id"] = ej.index.str.split("_").str[0].str.replace(":","_")+"_b38"
ej["gene"] = ej.index.str.split("_").str[-1].str.split("\\.").str[0]
ej.set_index(["variant_id", "gene"], inplace=True)
ej = ej.join(gpip[["pp_susie", "pp_fm", "ems_bin"]], how="left")
ejold["variant_id"] = ejold.variant_id.str.split(":").apply(lambda x: "_".join(x[:4]))+"_b38"
ejold["gene"] = ejold.gene_id.str.split("\\.").str[0]
ejold.set_index(["variant_id", "gene"], inplace=True)
ejold = ejold.join(gpip[["pp_susie", "pp_fm", "ems_bin"]], how="left")
#filter to defined ones:
ej = ej[~ej.ems_bin.isna()]
ejold = ejold[~ejold.ems_bin.isna()]
#save
ej.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eQTL_n1019_pip0001_gtexadded.tsv.gz", sep="\t", compression="gzip")
ejold.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eQTL_n465_pip0001_gtexadded.tsv.gz", sep="\t", compression="gzip")
#groupby:
emstb = ej.groupby(["min_pip_bin","ems_bin"]).size().unstack().fillna(0).astype(int)
emstb_old = ejold.groupby(["min_pip_bin","ems_bin"]).size().unstack().fillna(0).astype(int)
#gtex bin
ej["gtex_pip_bin"] = pd.Series(np.minimum(ej.pp_susie, ej.pp_fm)).apply(lambda x: bin_pip(x))
ejold["gtex_pip_bin"] = pd.Series(np.minimum(ejold.pp_susie, ejold.pp_fm)).apply(lambda x: bin_pip(x))
gtextb = ej.groupby(["min_pip_bin","gtex_pip_bin"]).size().unstack().fillna(0).astype(int)
gtextb_old = ejold.groupby(["min_pip_bin","gtex_pip_bin"]).size().unstack().fillna(0).astype(int)
#save these
emstb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eQTL_n1019_pip0001_emsbin.tsv", sep="\t")
emstb_old.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eQTL_n465_pip0001_emsbin.tsv", sep="\t")
gtextb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eQTL_n1019_pip0001_gtexbin.tsv", sep="\t")
gtextb_old.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eQTL_n465_pip0001_gtexbin.tsv", sep="\t")


#comparison with the previous version (n=1019 but old imputation panel)
#current version
e_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_finemap_pip0001.txt", sep=" ")
e_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1019_eqtl_susie_pip0001.txt", sep="\t")
e_fm["variant_id_hg38"] = e_fm.rsid.str.split("_").str[0]
e_fm["gene_id"] = e_fm.rsid.str.split("_").str[-1]
e_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
e_sus["variant_id_hg38"] = e_sus.rsid.str.split("_").str[0]
e_sus["gene_id"] = e_sus.rsid.str.split("_").str[-1]
e_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
ej = e_fm.join(e_sus.pip, how="outer")
ej.fillna(0, inplace=True) #anything lower than 000001 is 0
ej.columns = ["rsid", "pip_fm", "pip_sus"]
ej["fm_bin"] = ej.pip_fm.apply(lambda x: bin_pip(x))
ej["sus_bin"] = ej.pip_sus.apply(lambda x: bin_pip(x))
ej["min_pip"] = np.minimum(ej.pip_fm, ej.pip_sus)
ej["min_pip_bin"] = ej.min_pip.apply(lambda x: bin_pip(x))
#old version
e_fm_old = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/eqtl_finemap_pip0001.txt", sep=" ")
e_sus_old = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/eqtl_susie_pip0001.txt", sep="\t")
e_fm_old["variant_id_hg38"] = e_fm_old.rsid.str.split("_").str[0]
e_fm_old["gene_id"] = e_fm_old.rsid.str.split("_").str[-1]
e_fm_old.set_index(["variant_id_hg38","gene_id"], inplace=True)
e_sus_old["variant_id_hg38"] = e_sus_old.rsid.str.split("_").str[0]
e_sus_old["gene_id"] = e_sus_old.rsid.str.split("_").str[-1]
e_sus_old.set_index(["variant_id_hg38","gene_id"], inplace=True)
ej_old = e_fm_old.join(e_sus_old.pip, how="outer")
ej_old.fillna(0, inplace=True) #anything lower than 000001 is 0
ej_old.columns = ["rsid", "pip_fm", "pip_sus"]
ej_old["fm_bin"] = ej_old.pip_fm.apply(lambda x: bin_pip(x))
ej_old["sus_bin"] = ej_old.pip_sus.apply(lambda x: bin_pip(x))
ej_old["min_pip"] = np.minimum(ej_old.pip_fm, ej_old.pip_sus)
ej_old["min_pip_bin"] = ej_old.min_pip.apply(lambda x: bin_pip(x))

#join both:
ej = ej.join(ej_old, how="outer", rsuffix="_old")
ej.fillna(0, inplace=True)
tb = ej.groupby(["min_pip_bin_old", "min_pip_bin"]).size().unstack().fillna(0).astype(int)
tb.index = ["<0.001 or NA", "[0.001,0.01)", "[0.01,0.1)", "[0.1,0.5)", "[0.5,0.9)", "[0.9,1)", "==1"]
tb.columns = tb.index
tb.index.name = "PIP_old"
tb.columns.name = "PIP_updated"
#and add the total num variants:
e_fm_low = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/eqtl_finemap_n_below_pip0001.txt", sep='\s+', index_col=0)
e_sus_low = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/eqtl_susie_n_below_pip0001.txt", sep='\s+', index_col=0)#some random mistake and \t is not working as a tab...
tb.iloc[0,0] = e_fm_low.sum() - tb.sum(axis=1)[0] #all the fm bin==0 - (fm bin==0 and sus bin!=0) = (fm_bin==0 and sus_bin==0)
#and save for plotting
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/n1019_eqtl_pip_bef_vs_after_imputation.tsv", sep='\t')
#Also functional validation that current version is better
ej = ej.loc[:,["min_pip_bin","min_pip_bin_old"]]
ej["variant_id"] = ej.index.get_level_values(0).str.replace(":","_")+"_b38"
ej["gene"] = ej.index.get_level_values(1).str.split("\\.").str[0]
ej.set_index(["variant_id", "gene"], inplace=True)
#GTEx PIP and EMS
gpip = pd.read_csv("/Users/qingbowang/Desktop/resources/Whole_Blood_allpairs_gtex_pip.tsv.gz", sep="\t") #takes some time, large file
gpip["gene"] = gpip.gene_id.str.split("\\.").str[0] #takes some time..
gpip.set_index(["variant_id", "gene"], inplace=True)
ej = ej.join(gpip[["pp_susie", "pp_fm", "ems_bin"]], how="left")
#save
ej.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eQTL_comparison_imputation_gtexadded.tsv.gz", sep="\t", compression="gzip")




#plotting part
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

###p-value comparison
#vs n465:
plt.rcParams.update({'font.size': 12.5})
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019_vs_n465_stats.tsv", sep='\t', index_col=0)
tk = ["<2", "[2,3)", "[3,4)", "[4,5)", "[5,6)", "[6,7)", "[7,8)", "[8,9)", "[9,10)", "[10,20)", "[20,30)", "[30,50)", "[50,100)", "100$\leq$"]
tb.index = tk
tb.columns = tk
tb = tb.iloc[::-1,:] #to match the order
log_norm = LogNorm()
plt.figure(figsize=(8,8))
sns.heatmap(tb+1, annot=tb, fmt="d", square=True, linewidths=.5, norm=log_norm,
            cmap="viridis", cbar_kws={'label': 'count',
                                      "shrink": .6})
plt.yticks(rotation=35, fontsize=13)
plt.xticks(rotation=35, fontsize=13)
plt.xlabel("-log$_{10}$(P), previous version (n=465)", fontsize=15)
plt.ylabel("-log$_{10}$(P), this version (n=1019)", fontsize=15)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/vs_n465_eGenestats.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/vs_n465_eGenestats.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#vs n1019, pre-updating imputation panel:
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/minp_eqtl_n1019_vs_preimputation_stats.tsv", sep='\t', index_col=0)
tk = ["<2", "[2,3)", "[3,4)", "[4,5)", "[5,6)", "[6,7)", "[7,8)", "[8,9)", "[9,10)", "[10,20)", "[20,30)", "[30,50)", "[50,100)", "100$\leq$"]
tb.index = tk
tb.columns = tk
tb = tb.iloc[::-1,:] #to match the order
log_norm = LogNorm()
plt.figure(figsize=(8,8))
sns.heatmap(tb+1, annot=tb, fmt="d", square=True, linewidths=.5, norm=log_norm,
            cmap="viridis", cbar_kws={'label': 'count',
                                      "shrink": .6})
plt.yticks(rotation=35, fontsize=13)
plt.xticks(rotation=35, fontsize=13)
plt.xlabel("-log$_{10}$(P), previous imputation panel (n=1019)", fontsize=15)
plt.ylabel("-log$_{10}$(P), this version (n=1019)", fontsize=15)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/vs_n1019_prev_imputation_eGenestats.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/vs_n1019_prev_imputation_eGenestats.pdf', bbox_inches='tight', dpi=500)
plt.clf()

###PIP comparison
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/vs_n465_pip_updated.tsv", sep='\t', index_col=0)
c7 = ["lightskyblue", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
c9 = ["lightgray", "dimgray", "lightskyblue", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c9[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c9[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in previous version (n=465):", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pip_vs_n465_bar_full_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pip_vs_n465_bar_full_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#removing undefined and missing ones
tb = tb.iloc[:,2:]
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
#fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(9, 4.5), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
#ax[0].set_ylabel("N(variant-gene)", fontsize=15)
ax[0].set_ylabel("Count", fontsize=15)
ax[0].set_yticks([10**3,10**5,10**7])
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c7[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c7[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in previous version (n=465):", fontsize=12, title_fontsize=12)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pip_vs_n465_bar_main_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pip_vs_n465_bar_main_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#and the simple numbers increase:
c6 = ["tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
st = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/vs_n465_pip_number.tsv", sep='\t', index_col=0)
st.columns = ["This version (n=1,019)", "Previous version (n=465)"]
st.iloc[:,::-1].plot.bar(logy=True, color=["tab:gray", "tab:orange"], figsize=(5.2,2.8))
plt.xlabel("PIP bin")
plt.ylabel("Number of variant-genes")
plt.legend(title="Call set:")
tick_labels = plt.gca().get_xticklabels()
for i in range(len(tick_labels)):
    tick_labels[i].set_color(c6[i])
    tick_labels[i].set_rotation(30)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pip_vs_n465_numbers.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pip_vs_n465_numbers.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Functional evidence that updated version is better:
hts = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_updated.tsv", sep='\t', index_col=0)
hts["overall"] = hts.sum(axis=1)
htfrac = hts/hts.sum(axis=0)
htenr = (htfrac.T / htfrac["overall"]).T
hterr = (np.sqrt(htfrac*(1-htfrac) / hts.sum(axis=0)).T / htfrac.iloc[:,0]).T#分母のerrorは無視、as we did in EMS
#
hts = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_prev_missingvars_updated.tsv", sep='\t', index_col=0)
hts["overall"] = hts.sum(axis=1)
htfrac = hts/hts.sum(axis=0)
htenr2 = (htfrac.T / htfrac["overall"]).T
hterr2 = (np.sqrt(htfrac*(1-htfrac) / hts.sum(axis=0)).T / htfrac.iloc[:,0]).T#分母のerrorは無視、as we did in EMS
#
hts = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/tss_distance_distribution_prev_nonegenes_updated.tsv", sep='\t', index_col=0)
hts["overall"] = hts.sum(axis=1)
htfrac = hts/hts.sum(axis=0)
htenr3 = (htfrac.T / htfrac["overall"]).T
hterr3 = (np.sqrt(htfrac*(1-htfrac) / hts.sum(axis=0)).T / htfrac.iloc[:,0]).T#分母のerrorは無視、as we did in EMS
enr1 = htenr.applymap(lambda x: np.log2(x))
enr2 = htenr2.applymap(lambda x: np.log2(x))
enr3 = htenr3.applymap(lambda x: np.log2(x))
lower1 = (htenr-hterr).applymap(lambda x: np.log2(x))
lower2 = (htenr2-hterr2).applymap(lambda x: np.log2(x))
lower3 = (htenr3-hterr3).applymap(lambda x: np.log2(x))
upper1 = (htenr+hterr).applymap(lambda x: np.log2(x))
upper2 = (htenr2+hterr2).applymap(lambda x: np.log2(x))
upper3 = (htenr3+hterr3).applymap(lambda x: np.log2(x))

pips = ["PIP<0.001", "0.001$\leq$PIP<0.01", "0.01$\leq$PIP<0.1",
        "0.1$\leq$PIP<0.5", "0.5$\leq$PIP<0.9", "0.9$\leq$PIP<1", "PIP=1"]
fig, ax = plt.subplots(7, 1, sharex=True, figsize=(6, 6))
for i in range(7):
    ax[i].spines[['right', 'top']].set_visible(False)
    ax[i].axhline(y=0, linewidth=0.5, color="tab:gray", linestyle="--", zorder=-2)
    ax[i].errorbar(enr1.index-0.1, enr1.iloc[:,-2-i], [enr1.iloc[:,-2-i]-lower1.iloc[:,-2-i], upper1.iloc[:,-2-i]-enr1.iloc[:,-2-i]],
                   color=c7[-1-i], fmt="o", mec='black', mew=0.25, markersize=4)
    ax[i].errorbar(enr2.index, enr2.iloc[:, -2 - i],[enr2.iloc[:, -2 - i] - lower2.iloc[:, -2 - i], upper2.iloc[:, -2 - i] - enr2.iloc[:, -2 - i]],
                   color=c7[-1 - i], fmt="D", mec='black', mew=0.25, markersize=4)
    ax[i].errorbar(enr3.index+0.1, enr3.iloc[:, -2 - i],[enr3.iloc[:, -2 - i] - lower3.iloc[:, -2 - i], upper3.iloc[:, -2 - i] - enr3.iloc[:, -2 - i]],
                   color=c7[-1 - i], fmt="^", mec='black', mew=0.25, markersize=4)
for i in range(6):
    ax[i].text(6,enr1.iloc[:,-2-i].max()*0.9, pips[-1-i], color=c7[-1-i], ha="right")
ax[6].text(6,-0.22, pips[-1-6], color=c7[-1-6], ha="right")
ax[6].set_xlabel("Distance to TSS")
ax[3].set_ylabel("Log$_2$(enrichment)")
ax[6].set_xticks([2,3,4,5,6])
ax[6].set_xticklabels(["[$0,100$]","($10^{2.5},10^3$]","($10^{3.5},10^4$]","($10^{4.5},10^5$]","($10^{5.5},10^6$]"], rotation=30, fontsize=10)
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/tss_dist_vs_pip_updated.png", bbox_inches='tight', dpi=500)
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/tss_dist_vs_pip_updated.pdf", bbox_inches='tight', dpi=500)
plt.clf()

#functional evidence in terms of EMS and GTEx:
emstb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eQTL_n1019_pip0001_emsbin.tsv", sep="\t", index_col=0)
emstb_old = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eQTL_n465_pip0001_emsbin.tsv", sep="\t", index_col=0)
gtextb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eQTL_n1019_pip0001_gtexbin.tsv", sep="\t", index_col=0)
gtextb_old = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eQTL_n465_pip0001_gtexbin.tsv", sep="\t", index_col=0)

emstb_norm = (emstb.T/emstb.sum(axis=1)).T
emstb_old_norm = (emstb_old.T/emstb_old.sum(axis=1)).T
gtextb_norm = (gtextb.T/gtextb.sum(axis=1)).T
gtextb_old_norm = (gtextb_old.T/gtextb_old.sum(axis=1)).T
#EMS:
tk = ["[0.001,0.01)","[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
emstb_norm.columns = ["(0.01,0.1]","(0.1,1]","(1,10]","(10,100]","(100,1000]","1000<"]
from matplotlib import cm
vir = cm.viridis
fig, ax = plt.subplots(figsize=(6,3))
handles = []#legend handles for EMS
#fill the edge first
ax0 = ax.bar(np.arange(emstb_norm.shape[0])-0.21, [1]*emstb_norm.shape[0], fill = False, edgecolor = "tab:gray", linewidth=1.5, width=0.4, label="Previous version\n(n=465; left)")
ax1 = ax.bar(np.arange(emstb_norm.shape[0])+0.21, [1]*emstb_norm.shape[0], fill = False, edgecolor = "navy", linewidth=1.5, width=0.4, label="This version\n(n=1,019; right)")
first_legend = ax.legend(handles=[ax0, ax1], bbox_to_anchor=(1.01,0.9), loc='center left', fontsize=11, title_fontsize=11)
ax.add_artist(first_legend)
axn = ax.bar(np.arange(emstb_norm.shape[0])+0.2+0.01, emstb_norm.iloc[:,0], color=vir(int(0/(emstb_norm.shape[1]-1)*256)), edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005,label=emstb_norm.columns[0])
handles.append(axn)
for i in range(1,emstb_norm.shape[1]):
        axn = ax.bar(np.arange(emstb_norm.shape[0])+0.2+0.01, emstb_norm.iloc[:,i], bottom=emstb_norm.iloc[:, :i].sum(axis=1),label=emstb_norm.columns[i],
                      color=vir(int(i/(emstb_norm.shape[1]-1)*256)),edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
        handles.append(axn)
ax.bar(np.arange(emstb_old_norm.shape[0])-0.2-0.01, emstb_old_norm.iloc[:,0], color=vir(int(0/(emstb_old_norm.shape[1]-1)*256)), edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
for i in range(1,emstb_old_norm.shape[1]):
        ax.bar(np.arange(emstb_old_norm.shape[0])-0.2-0.01, emstb_old_norm.iloc[:,i], bottom=emstb_old_norm.iloc[:, :i].sum(axis=1),
                      color=vir(int(i/(emstb_old_norm.shape[1]-1)*256)),edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
ax.legend(handles=handles, bbox_to_anchor=(1.01,0.3), loc='center left', title="Expression modifier\nscore (EMS)", fontsize=11, title_fontsize=11)
ax.set_xlabel("PIP bin")
ax.set_ylabel("Fraction")
ax.set_xlim([-1+0.5,emstb_norm.shape[0]-0.5])
ax.set_ylim([0,1])
ax.set_xticks(np.arange(emstb_norm.shape[0]), tk, rotation=30)
tick_labels = plt.gca().get_xticklabels()
for i in range(len(tick_labels)):
    tick_labels[i].set_color(c6[i])
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/n1019_vs_n465_ems.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/plots/n1019_vs_n465_ems.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#GTEx:
c7 = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tk = ["[0.001,0.01)","[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
gtextb_norm.columns = ["<0.001", "[0.001,0.01)", "[0.01,0.1)", "[0.1,0.5)", "[0.5,0.9)", "[0.9,1)", "1"]
fig, ax = plt.subplots(figsize=(6,3))
handles = []#legend handles for EMS
#fill the edge first
ax0 = ax.bar(np.arange(gtextb_norm.shape[0])-0.21, [1]*gtextb_norm.shape[0], fill = False, edgecolor = "tab:gray", linewidth=1.5, width=0.4, label="Previous version \n(n=465; left)")
ax1 = ax.bar(np.arange(gtextb_norm.shape[0])+0.21, [1]*gtextb_norm.shape[0], fill = False, edgecolor = "navy", linewidth=1.5, width=0.4, label="This version \n(n=1,019; right)")
first_legend = ax.legend(handles=[ax0, ax1], bbox_to_anchor=(1.01,0.9), loc='center left', fontsize=11, title_fontsize=11)
ax.add_artist(first_legend)
axn = ax.bar(np.arange(gtextb_norm.shape[0])+0.2+0.01, gtextb_norm.iloc[:,0], color=c7[0], edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005,label=gtextb_norm.columns[0])
handles.append(axn)
for i in range(1,gtextb_norm.shape[1]):
        axn = ax.bar(np.arange(gtextb_norm.shape[0])+0.2+0.01, gtextb_norm.iloc[:,i], bottom=gtextb_norm.iloc[:, :i].sum(axis=1),label=gtextb_norm.columns[i],
                      color=c7[i],edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
        handles.append(axn)
ax.bar(np.arange(gtextb_old_norm.shape[0])-0.2-0.01, gtextb_old_norm.iloc[:,0], color=c7[0], edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
for i in range(1,gtextb_old_norm.shape[1]):
        ax.bar(np.arange(gtextb_old_norm.shape[0])-0.2-0.01, gtextb_old_norm.iloc[:,i], bottom=gtextb_old_norm.iloc[:, :i].sum(axis=1),
                      color=c7[i],edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
ax.legend(handles=handles, bbox_to_anchor=(1.01,0.3), loc='center left', title="PIP in GTEx", fontsize=11, title_fontsize=11)
ax.set_xlabel("PIP bin")
ax.set_ylabel("Fraction")
ax.set_xlim([-1+0.5,gtextb_norm.shape[0]-0.5])
ax.set_ylim([0,1])
ax.set_xticks(np.arange(gtextb_norm.shape[0]), tk, rotation=30)
tick_labels = plt.gca().get_xticklabels()
for i in range(len(tick_labels)):
    tick_labels[i].set_color(c6[i])
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/n1019_vs_n465_gtex.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/plots/n1019_vs_n465_gtex.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#SuSiE vs FM
tk = ["[0,0.001)","[0.001,0.01)","[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/eqtl_sus_vs_fm_bins.tsv", sep='\t', index_col=0)
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
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pip_fm_vs_susie.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pip_fm_vs_susie.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#before and after imputation:
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/n1019_eqtl_pip_bef_vs_after_imputation.tsv", sep='\t', index_col=0)
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
plt.title("eQTL PIP before and after updating imputation panel")
plt.xlabel("PIP_old", fontsize=16)
plt.ylabel("PIP_updated", fontsize=16)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/n1019_eqtl_pip_bef_vs_after_imputation_heatmap.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/n1019_eqtl_pip_bef_vs_after_imputation_heatmap.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#and fraction:
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/fmoutputs/summary/n1019_eqtl_pip_bef_vs_after_imputation.tsv", sep='\t', index_col=0)
tb.index = ["[0,0.001)\n or NA","[0.001,0.01)","[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
tb.columns = tb.index
tbnorm = (tb.T/tb.sum(axis=1)).T
c7 = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tbnorm.plot.bar(stacked=True, color=c7, figsize=(7,3.5)).legend(title="PIP_updated", loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.title("eQTL PIP after updating imputation panel")
plt.xlabel("PIP_old")
plt.ylabel("Fraction")
plt.xticks(rotation=30)
for i in range(7):
    plt.gca().get_xticklabels()[i].set_color(c7[i])
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/pip_updated_by_imputation_eqtl.png", dpi=400)
plt.savefig("/Users/qingbowang/Desktop/plots/pip_updated_by_imputation_eqtl.pdf", dpi=400)
plt.clf()

#and GTEx enrichment:
ej = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eQTL_comparison_imputation_gtexadded.tsv.gz",
                 sep="\t", compression="gzip", index_col=[0,1])
ex = ej[~ej.ems_bin.isna()]
def bin_pip_rough2(pip_bin):
    #if pip_bin >= 0.5:return 1  # more than 0.5 => "high"
    if pip_bin>=0.9: return 1 #more than 0.9 => "high"
    elif pip_bin<0.5: return 0 #less than 0.5 => "low"
    #elif pip_bin < 0.1:return 0  # less than 0.5 => "low"
    else: return (np.nan)
ex["pip_old_bin2"] = ex.min_pip_bin_old.apply(lambda x: bin_pip_rough2(x))
ex["pip_new_bin2"] = ex.min_pip_bin.apply(lambda x: bin_pip_rough2(x))
ex["pip_min_gtex"] = np.minimum(ex.pp_fm, ex.pp_susie)
ex["fm_gtex_bin2"] = ex.pip_min_gtex.apply(lambda x: bin_pip(x))
comparison = ex.groupby(["pip_old_bin2", "pip_new_bin2", "fm_gtex_bin2"]).size().unstack()
comparison = comparison.loc[~comparison.index.get_level_values(0).isna(),:]
comparison = comparison.loc[~comparison.index.get_level_values(1).isna(),:]
comparison = comparison.loc[abs(comparison.index.get_level_values(0)-comparison.index.get_level_values(1))==1,:]
comparison_norm = (comparison.T/comparison.sum(axis=1)).T #わかりにくいか..
comparison2 = ex.groupby(["pip_old_bin2", "pip_new_bin2", "ems_bin"]).size().unstack()
comparison2 = comparison2.loc[comparison2.index.get_level_values(0)!=100,:]
comparison2 = comparison2.loc[comparison2.index.get_level_values(1)!=100,:]
comparison2 = comparison2.loc[abs(comparison2.index.get_level_values(0)-comparison2.index.get_level_values(1))==1,:]
comparison2_norm = (comparison2.T/comparison2.sum(axis=1)).T #わかりにくいか..
print (comparison)
print (comparison_norm)
print (comparison2)
print (comparison2_norm)
#plot these:
from matplotlib import cm
vir = cm.viridis
virs = []
for i in range(6):
    virs.append(vir(int(i/5 * 256)))
c7 = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
comparison_norm.columns = ["<0.001", "[0.001,0.01)", "[0.01,0.1)", "[0.1,0.5)", "[0.5,0.9)", "[0.9,1)", "1"]
comparison_norm.plot.bar(stacked=True, color=c7, figsize=(4,4), width=.9).legend(title="PIP in GTEx:", loc='center left',bbox_to_anchor=(1.0, 0.5), fontsize=11, title_fontsize=11)
plt.xlabel("PIP difference direction")
plt.ylabel("Fraction")
plt.ylim([0,1])
plt.xticks([0,1],["From: Low\nTo: High", "From: High\nTo: Low"], rotation=45)
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/pip_updated_by_imputation_eqtl_highlight.png", dpi=400)
plt.savefig("/Users/qingbowang/Desktop/plots/pip_updated_by_imputation_eqtl_highlight.pdf", dpi=400)
plt.clf()

comparison2_norm.columns = ["(0.01,0.1]", "(0.1,1]","(1,10]", "(10,100]", "(100,1000]", "1000<"]
comparison2_norm.plot.bar(stacked=True, color=virs, figsize=(4,4), width=.9).legend(title="EMS", loc='center left',bbox_to_anchor=(1.0, 0.5), fontsize=11, title_fontsize=11)
plt.xlabel("PIP difference direction")
plt.tight_layout()
plt.ylabel("Fraction")
plt.ylim([0,1])
plt.xticks([0,1],["From: Low\nTo: High", "From: High\nTo: Low"], rotation=45)
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/pip_updated_by_imputation_eqtl_highlight2.png", dpi=400)
plt.savefig("/Users/qingbowang/Desktop/plots/pip_updated_by_imputation_eqtl_highlight2.pdf", dpi=400)
plt.clf()

