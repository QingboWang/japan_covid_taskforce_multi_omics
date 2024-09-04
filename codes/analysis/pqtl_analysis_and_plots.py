import pandas as pd
import numpy as np
import time as tm

#min(p) per gene
minp = []
for chr in list(range(1,22+1))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(chr)
    df = pd.read_csv(input_file, sep="\t")
    minp.append(df.groupby("gene_id").pval_nominal.min())
    print ("done {0}, {1}".format(chr, tm.ctime()))
minp = pd.concat(minp)
minp[minp==0] = minp[minp>0].min()
minp.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/pqtl_minp_pergene.txt", sep="\t")

#SuSiE vs FINEMAP
p_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_finemap_pip0001.txt", sep=' ')
p_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_susie_pip0001.txt", sep='\t')
p_fm["variant_id_hg38"] = p_fm.rsid.str.split("_").str[0]
p_fm["gene_id"] = p_fm.rsid.str.split("_").str[-1]
p_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
p_sus["variant_id_hg38"] = p_sus.rsid.str.split("_").str[0]
p_sus["gene_id"] = p_sus.rsid.str.split("_").str[-1]
p_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
pj = p_fm.join(p_sus.pip, how="outer")
pj.fillna(0, inplace=True) #anything lower than 0001 is 0
pj.columns = ["rsid", "pip_fm", "pip_sus"]
def bin_pip(pip):
    if pip<0.001: return 0
    elif pip<0.01: return 0.001
    elif pip<0.1: return 0.01
    elif pip<0.5: return 0.1
    elif pip<0.9: return 0.5
    elif pip<1: return 0.9
    elif pip==1: return 1
    else: return np.nan
pj["fm_bin"] = pj.pip_fm.apply(lambda x: bin_pip(x))
pj["sus_bin"] = pj.pip_sus.apply(lambda x: bin_pip(x))
pj["min_pip"] = np.minimum(pj.pip_fm, pj.pip_sus)
pj["min_pip_bin"] = pj.min_pip.apply(lambda x: bin_pip(x))
#table:
tb = pj.groupby(["fm_bin", "sus_bin"]).size().unstack().fillna(0).astype(int) #very good agreement confirmed (so far)
#also add number for 0,0:
p_fm_low = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_finemap_n_below_pip0001.txt", sep='\t', index_col=0)
p_sus_low = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_susie_n_below_pip0001.txt", sep='\t', index_col=0)
tb.iloc[0,0] = p_fm_low.sum() - tb.sum(axis=1)[0] #subtracting those where fm bin != 0
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/pqtl_sus_vs_fm_bins.tsv", sep='\t')



##Comparison with Zhang et al
#Reading their data and matching the ensgIDs to gencode v30
ensgids = pd.read_csv("~/Desktop/taskforce_n1102/gencode.v30.genes.parsed.tsv", sep='\t')
ensgids.index =ensgids.gene_name
aa = pd.read_csv("~/Desktop/resources/pqtl_aa_all.tsv.gz", sep="\t")
aa.index = aa.gene_name
aa = aa.join(ensgids.ensg_id, how="left")
ea = pd.read_csv("~/Desktop/resources/pqtl_ea_all.tsv.gz", sep="\t")
ea.index = ea.gene_name
ea = ea.join(ensgids.ensg_id, how="left")
#Set the variant ID properly:
aa["variant_hg38"] = "chr"+aa.iloc[:,0].astype(str)+":"+aa.iloc[:,1].astype(str)+":"+aa.iloc[:,3]+":"+aa.iloc[:,4]
ea["variant_hg38"] = "chr"+ea.iloc[:,0].astype(str)+":"+ea.iloc[:,1].astype(str)+":"+ea.iloc[:,3]+":"+ea.iloc[:,4]
#QC: remove ensg_id==NAs here
sum(aa.ensg_id.isna())/aa.shape[0] #2% is not low... but seems fine.
aa[aa.ensg_id.isna()]
sum(ea.ensg_id.isna())/ea.shape[0] #2% is not low... but seems fine. Just remove them (although some of them are obvious alias..)
aa = aa[~aa.ensg_id.isna()]
ea = ea[~ea.ensg_id.isna()]

#Gene level tested in JCTF flag:
aa.index = aa.ensg_id
ea.index = ea.ensg_id
genes_p = pd.read_csv("~/Desktop/taskforce_n1102/n1369.protein.expression.bed", sep='\t', index_col=3)
#genes_p = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/protein_measured_genelist.tsv", sep='\t', header=None, index_col=0)#old
genes_p["is_tested_protein"] = True
aa = aa.join(genes_p.is_tested_protein, how="left")
ea = ea.join(genes_p.is_tested_protein, how="left")
aa["is_tested_protein"].fillna(False, inplace=True)
ea["is_tested_protein"].fillna(False, inplace=True)

#Variant level "tested in JCTF" flag:
vars = []
for chr in list(range(1,22+1))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(chr)
    df = pd.read_csv(input_file, sep='\t')
    vars = vars + list(df.variant_id.str.split("_").str[0].unique())
    print ("done chr{0}, {1}".format(chr, tm.ctime()))
vars = pd.Series(vars)
vars.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/n1384_protein_measured_variantlist.tsv", sep='\t', index=False)
vars = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/n1384_protein_measured_variantlist.tsv", sep='\t', header=None, index_col=0)
vars["is_tested_variant"] = True
#add that flag:
aa.index = aa.variant_hg38
ea.index = ea.variant_hg38
aa = aa.join(vars, how="left")
ea = ea.join(vars, how="left")
aa["is_tested_variant"].fillna(False, inplace=True)
ea["is_tested_variant"].fillna(False, inplace=True)
#sanity check
aa.groupby(["is_tested_protein", "is_tested_variant"]).size().unstack()
#add PIP when >0.001
p_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_finemap_pip0001.txt", sep=' ')
p_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_susie_pip0001.txt", sep='\t')
#actually we only need PIP>=0.001? No we need till 0.001
p_fm = p_fm[p_fm.prob>=0.001]
p_fm["variant_id_hg38"] = p_fm.rsid.str.split("_").str[0]
p_fm["gene_id"] = p_fm.rsid.str.split("_").str[-1]
p_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
p_sus = p_sus[p_sus.pip>=0.001]
p_sus["variant_id_hg38"] = p_sus.rsid.str.split("_").str[0]
p_sus["gene_id"] = p_sus.rsid.str.split("_").str[-1]
p_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
pj = p_fm.join(p_sus.pip, how="outer")
pj.fillna(0, inplace=True) #anything lower than 0.001 is 0
pj.columns = ["rsid", "p_pip_fm", "p_pip_sus"]
pj["p_min_pip"] = np.minimum(pj.p_pip_fm, pj.p_pip_sus)
aa.set_index(["variant_hg38","ensg_id"], inplace=True)
ea.set_index(["variant_hg38","ensg_id"], inplace=True)
pj.index.names = ["variant_hg38","ensg_id"]
aa = aa.join(pj[["p_pip_fm", "p_pip_sus","p_min_pip"]], how="left")
ea = ea.join(pj[["p_pip_fm", "p_pip_sus","p_min_pip"]], how="left")
#Now compare the two.
#first their aa vs our SuSiE PIP
aa["susie_pip_bin"] = 0
aa.loc[aa.p_pip_sus>0.001, "susie_pip_bin"] = 0.001
aa.loc[aa.p_pip_sus>0.01, "susie_pip_bin"] = 0.01
aa.loc[aa.p_pip_sus>0.1, "susie_pip_bin"] = 0.1
aa.loc[aa.p_pip_sus>0.5, "susie_pip_bin"] = 0.5
aa.loc[aa.p_pip_sus>0.9, "susie_pip_bin"] = 0.9
aa.loc[aa.p_pip_sus==1, "susie_pip_bin"] = 1
aa["missing_status"] = "Exists"
aa.loc[(~aa.is_tested_protein) & (~aa.is_tested_variant),"missing_status"] = "Both missing"
aa.loc[(~aa.is_tested_protein) & (aa.is_tested_variant),"missing_status"] = "Protein missing"
aa.loc[(aa.is_tested_protein) & (~aa.is_tested_variant),"missing_status"] = "Variant missing"
aafilt = aa[aa.is_tested_protein & aa.is_tested_variant]
tb = []
tbmis = []
for i in np.arange(101) / 100:
    j = aafilt[aafilt.iloc[:,14] >= i].susie_pip_bin.value_counts()
    j2 = aa[aa.iloc[:, 14] >= i].missing_status.value_counts()
    tb.append(j)
    tbmis.append(j2)
    print ("done {0}, {1}".format(i, tm.ctime()))
tb = pd.DataFrame(tb)
tb.index = range(101)
tb.columns = ["(0,0.001]","(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tbmis = pd.DataFrame(tbmis)
tbmis.index = range(101)
del tbmis["Exists"]
tb = pd.concat([tb, tbmis], axis=1).fillna(0).astype(int)
print(tb)
print(tb.apply(lambda x: x / sum(x), axis=1))
print(tb.iloc[:,:-3].apply(lambda x: x / sum(x), axis=1)) #filtering to existing ones
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/lit_aa_vs_our_pqtl_susie_stats.tsv", sep='\t')
#we are much more conservative
#the min(sus, fm) we are using (probably will be even more conservative):
aa["min_pip_bin"] = 0
aa.loc[aa.p_min_pip>0.001, "min_pip_bin"] = 0.001
aa.loc[aa.p_min_pip>0.01, "min_pip_bin"] = 0.01
aa.loc[aa.p_min_pip>0.1, "min_pip_bin"] = 0.1
aa.loc[aa.p_min_pip>0.5, "min_pip_bin"] = 0.5
aa.loc[aa.p_min_pip>0.9, "min_pip_bin"] = 0.9
aa.loc[aa.p_min_pip==1, "min_pip_bin"] = 1
aafilt = aa[aa.is_tested_protein & aa.is_tested_variant]
tb = []
tbmis = []
for i in np.arange(101) / 100:
    j = aafilt[aafilt.iloc[:,14] >= i].min_pip_bin.value_counts()
    j2 = aa[aa.iloc[:, 14] >= i].missing_status.value_counts() #redundant but doing again just in case
    tb.append(j)
    tbmis.append(j2)
    print ("done {0}, {1}".format(i, tm.ctime()))
tb = pd.DataFrame(tb)
tb.index = range(101)
tb.columns = ["(0,0.001]","(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tbmis = pd.DataFrame(tbmis)
tbmis.index = range(101)
del tbmis["Exists"]
tb = pd.concat([tb, tbmis], axis=1).fillna(0).astype(int)
print(tb)
print(tb.apply(lambda x: x / sum(x), axis=1))
print(tb.iloc[:,:-3].apply(lambda x: x / sum(x), axis=1)) #filtering to existing ones
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/lit_aa_vs_our_pqtl_minpip_stats.tsv", sep='\t')

#same for ea:
ea["susie_pip_bin"] = 0
ea.loc[ea.p_pip_sus>0.001, "susie_pip_bin"] = 0.001
ea.loc[ea.p_pip_sus>0.01, "susie_pip_bin"] = 0.01
ea.loc[ea.p_pip_sus>0.1, "susie_pip_bin"] = 0.1
ea.loc[ea.p_pip_sus>0.5, "susie_pip_bin"] = 0.5
ea.loc[ea.p_pip_sus>0.9, "susie_pip_bin"] = 0.9
ea.loc[ea.p_pip_sus==1, "susie_pip_bin"] = 1

ea["missing_status"] = "Exists"
ea.loc[(~ea.is_tested_protein) & (~ea.is_tested_variant),"missing_status"] = "Both missing"
ea.loc[(~ea.is_tested_protein) & (ea.is_tested_variant),"missing_status"] = "Protein missing"
ea.loc[(ea.is_tested_protein) & (~ea.is_tested_variant),"missing_status"] = "Variant missing"

eafilt = ea[ea.is_tested_protein & ea.is_tested_variant]
tb = []
tbmis = []
for i in np.arange(101) / 100:
    j = eafilt[eafilt.iloc[:,14] >= i].susie_pip_bin.value_counts()
    j2 = ea[ea.iloc[:, 14] >= i].missing_status.value_counts()
    tb.append(j)
    tbmis.append(j2)
    print ("done {0}, {1}".format(i, tm.ctime()))
tb = pd.DataFrame(tb)
tb.index = range(101)
tb.columns = ["(0,0.001]","(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tbmis = pd.DataFrame(tbmis)
tbmis.index = range(101)
del tbmis["Exists"]
tb = pd.concat([tb, tbmis], axis=1).fillna(0).astype(int)
print(tb)
print(tb.apply(lambda x: x / sum(x), axis=1))
print(tb.iloc[:,:-3].apply(lambda x: x / sum(x), axis=1)) #filtering to existing ones
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/lit_ea_vs_our_pqtl_susie_stats.tsv", sep='\t')
#we are much more conservative

#the min(sus, fm) we are using (probably will be even more conservative):
ea["min_pip_bin"] = 0
ea.loc[ea.p_min_pip>0.001, "min_pip_bin"] = 0.001
ea.loc[ea.p_min_pip>0.01, "min_pip_bin"] = 0.01
ea.loc[ea.p_min_pip>0.1, "min_pip_bin"] = 0.1
ea.loc[ea.p_min_pip>0.5, "min_pip_bin"] = 0.5
ea.loc[ea.p_min_pip>0.9, "min_pip_bin"] = 0.9
ea.loc[ea.p_min_pip==1, "min_pip_bin"] = 1
eafilt = ea[ea.is_tested_protein & ea.is_tested_variant]
tb = []
tbmis = []
for i in np.arange(101) / 100:
    j = eafilt[eafilt.iloc[:,14] >= i].min_pip_bin.value_counts()
    j2 = ea[ea.iloc[:, 14] >= i].missing_status.value_counts() #redundant but doing again just in case
    tb.append(j)
    tbmis.append(j2)
    print ("done {0}, {1}".format(i, tm.ctime()))
tb = pd.DataFrame(tb)
tb.index = range(101)
tb.columns = ["(0,0.001]","(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tbmis = pd.DataFrame(tbmis)
tbmis.index = range(101)
del tbmis["Exists"]
tb = pd.concat([tb, tbmis], axis=1).fillna(0).astype(int)
print(tb)
print(tb.apply(lambda x: x / sum(x), axis=1))
print(tb.iloc[:,:-3].apply(lambda x: x / sum(x), axis=1)) #filtering to existing ones
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/lit_ea_vs_our_pqtl_minpip_stats.tsv", sep='\t')

#also would like to do the other way around:
#i.e. when exists, is our result replicated in at least one of two?
pj = pj.join(aa.iloc[:,14], how="left")
pj.rename(columns = {"14":"pip_aa"}, inplace=True)
pj.pip_aa.fillna(-1, inplace=True)#anything that is not defined let it be -1
pj = pj.join(ea.iloc[:,14], how="left")
pj.rename(columns = {"14":"pip_ea"}, inplace=True)
pj.pip_ea.fillna(-1, inplace=True)#anything that is not defined let it be -1

#also classify the missing; variant missing, gene missing, or both
aa_genes = aa.index.get_level_values(1).unique()
ea_genes = ea.index.get_level_values(1).unique()
aa_vars = aa.index.get_level_values(0).unique()
ea_vars = ea.index.get_level_values(0).unique()
either_genes = aa_genes.union(ea_genes)
either_vars = aa_vars.union(ea_vars)
aa_genes = pd.DataFrame({"aa_protein_tested":True}, index=aa_genes)
ea_genes = pd.DataFrame({"ea_protein_tested":True}, index=ea_genes)
either_genes = pd.DataFrame({"either_protein_tested":True}, index=either_genes)
aa_vars = pd.DataFrame({"aa_var_tested":True}, index=aa_vars)
ea_vars = pd.DataFrame({"ea_var_tested":True}, index=ea_vars)
either_vars = pd.DataFrame({"either_var_tested":True}, index=either_vars)
pj = pj.join(aa_genes, how="left")
pj = pj.join(ea_genes, how="left")
pj = pj.join(either_genes, how="left")
pj = pj.join(aa_vars, how="left")
pj = pj.join(ea_vars, how="left")
pj = pj.join(either_vars, how="left")
for col in pj.columns[-6:]:
    pj.loc[:,col].fillna(False, inplace=True)
pj["missingness_class_aa"] = "Exists"
pj.loc[~pj.aa_var_tested, "missingness_class_aa"] = "Variant missing"
pj.loc[~pj.aa_protein_tested, "missingness_class_aa"] = "Gene missing"
pj.loc[(~pj.aa_protein_tested)&(~pj.aa_var_tested), "missingness_class_aa"] = "Gene & Variant\nmissing"
pj["missingness_class_ea"] = "Exists"
pj.loc[~pj.ea_var_tested, "missingness_class_ea"] = "Variant missing"
pj.loc[~pj.ea_protein_tested, "missingness_class_ea"] = "Gene missing"
pj.loc[(~pj.ea_protein_tested)&(~pj.ea_var_tested), "missingness_class_ea"] = "Gene & Variant\nmissing"
pj["missingness_class_either"] = "Exists"
pj.loc[~pj.either_var_tested, "missingness_class_either"] = "Variant missing"
pj.loc[~pj.either_protein_tested, "missingness_class_either"] = "Gene missing"
pj.loc[(~pj.either_protein_tested)&(~pj.either_var_tested), "missingness_class_either"] = "Gene & Variant\nmissing"


pj["aa_pip_bin"] = 0 #-1
pj.loc[pj.pip_aa>0.001, "aa_pip_bin"] = 0.001
pj.loc[pj.pip_aa>0.01, "aa_pip_bin"] = 0.01
pj.loc[pj.pip_aa>0.1, "aa_pip_bin"] = 0.1
pj.loc[pj.pip_aa>0.5, "aa_pip_bin"] = 0.5
pj.loc[pj.pip_aa>0.9, "aa_pip_bin"] = 0.9
pj.loc[pj.pip_aa==1, "aa_pip_bin"] = 1

pj["ea_pip_bin"] = 0#-1
pj.loc[pj.pip_ea>0.001, "ea_pip_bin"] = 0.001
pj.loc[pj.pip_ea>0.01, "ea_pip_bin"] = 0.01
pj.loc[pj.pip_ea>0.1, "ea_pip_bin"] = 0.1
pj.loc[pj.pip_ea>0.5, "ea_pip_bin"] = 0.5
pj.loc[pj.pip_ea>0.9, "ea_pip_bin"] = 0.9
pj.loc[pj.pip_ea==1, "ea_pip_bin"] = 1

pj["aaea_max_pip_bin"] = 0#-1
pj.loc[(pj.pip_aa>0.001)|(pj.pip_ea>0.001), "aaea_max_pip_bin"] = 0.001
pj.loc[(pj.pip_aa>0.01)|(pj.pip_ea>0.01), "aaea_max_pip_bin"] = 0.01
pj.loc[(pj.pip_aa>0.1)|(pj.pip_ea>0.1), "aaea_max_pip_bin"] = 0.1
pj.loc[(pj.pip_aa>0.5)|(pj.pip_ea>0.5), "aaea_max_pip_bin"] = 0.5
pj.loc[(pj.pip_aa>0.9)|(pj.pip_ea>0.9), "aaea_max_pip_bin"] = 0.9
pj.loc[(pj.pip_aa==1)|(pj.pip_ea==1), "aaea_max_pip_bin"] = 1

pjexists_aa = pj[pj.aa_var_tested & pj.aa_protein_tested] #for aa
pjmissing_aa = pj[~(pj.aa_var_tested & pj.aa_protein_tested)] #for aa
pjexists_ea = pj[pj.ea_var_tested & pj.ea_protein_tested] #for ea
pjmissing_ea = pj[~(pj.ea_var_tested & pj.ea_protein_tested)] #for ea
pjexists_either = pj[pj.either_var_tested & pj.either_protein_tested] #for either
pjmissing_either = pj[~(pj.either_var_tested & pj.either_protein_tested)] #for either
tbaa = []
tbea = []
tbeither = []
tb2aa = []
tb2ea = []
tb2either = [] #2 is for the missing ones
for i in np.arange(101) / 100:
    sub = pjexists_aa[pjexists_aa.p_pip_sus >= i] #for aa
    j = sub.aa_pip_bin.value_counts()
    tbaa.append(j)
    sub = pjmissing_aa[pjmissing_aa.p_pip_sus >= i]
    j2 = sub.missingness_class_aa.value_counts()
    tb2aa.append(j2)
    sub = pjexists_ea[pjexists_ea.p_pip_sus >= i] #for ea
    j = sub.ea_pip_bin.value_counts()
    tbea.append(j)
    sub = pjmissing_ea[pjmissing_ea.p_pip_sus >= i]
    j2 = sub.missingness_class_ea.value_counts()
    tb2ea.append(j2)
    sub = pjexists_either[pjexists_either.p_pip_sus >= i] #for either
    j = sub.aaea_max_pip_bin.value_counts()
    tbeither.append(j)
    sub = pjmissing_either[pjmissing_either.p_pip_sus >= i]
    j2 = sub.missingness_class_either.value_counts()
    tb2either.append(j2)
    print ("done {0}, {1}".format(i, tm.ctime()))
tbaa = pd.DataFrame(tbaa).fillna(0).astype(int)
tbaa.index = range(101)
tbaa.columns = ["[0,0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tb2aa = pd.DataFrame(tb2aa).fillna(0).astype(int)
tb2aa.index = range(101)
tbaa = pd.concat([tb2aa, tbaa], axis=1)
tbaa.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_suspip_vs_lit_aa_stats.tsv", sep='\t')
tbea = pd.DataFrame(tbea).fillna(0).astype(int)
tbea.index = range(101)
tbea.columns = ["[0,0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tb2ea = pd.DataFrame(tb2ea).fillna(0).astype(int)
tb2ea.index = range(101)
tbea = pd.concat([tb2ea, tbea], axis=1)
tbea.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_suspip_vs_lit_ea_stats.tsv", sep='\t')
tbeither = pd.DataFrame(tbeither).fillna(0).astype(int)
tbeither.index = range(101)
tbeither.columns = ["[0,0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tb2either = pd.DataFrame(tb2either).fillna(0).astype(int)
tb2either.index = range(101)
tbeither = pd.concat([tb2either, tbeither], axis=1)
tbeither.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_suspip_vs_lit_either_stats.tsv", sep='\t')
#also versus min pip instead of SuSiE:
tbaa = []
tbea = []
tbeither = []
tb2aa = []
tb2ea = []
tb2either = [] #2 is for the missing ones
for i in np.arange(101) / 100:
    sub = pjexists_aa[pjexists_aa.p_min_pip >= i] #for aa
    j = sub.aa_pip_bin.value_counts()
    tbaa.append(j)
    sub = pjmissing_aa[pjmissing_aa.p_min_pip >= i]
    j2 = sub.missingness_class_aa.value_counts()
    tb2aa.append(j2)
    sub = pjexists_ea[pjexists_ea.p_min_pip >= i] #for ea
    j = sub.ea_pip_bin.value_counts()
    tbea.append(j)
    sub = pjmissing_ea[pjmissing_ea.p_min_pip >= i]
    j2 = sub.missingness_class_ea.value_counts()
    tb2ea.append(j2)
    sub = pjexists_either[pjexists_either.p_min_pip >= i] #for either
    j = sub.aaea_max_pip_bin.value_counts()
    tbeither.append(j)
    sub = pjmissing_either[pjmissing_either.p_min_pip >= i]
    j2 = sub.missingness_class_either.value_counts()
    tb2either.append(j2)
    print ("done {0}, {1}".format(i, tm.ctime()))
tbaa = pd.DataFrame(tbaa).fillna(0).astype(int)
tbaa.index = range(101)
tbaa.columns = ["[0,0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tb2aa = pd.DataFrame(tb2aa).fillna(0).astype(int)
tb2aa.index = range(101)
tbaa = pd.concat([tb2aa, tbaa], axis=1)
tbaa.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_minpip_vs_lit_aa_stats.tsv", sep='\t')
tbea = pd.DataFrame(tbea).fillna(0).astype(int)
tbea.index = range(101)
tbea.columns = ["[0,0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tb2ea = pd.DataFrame(tb2ea).fillna(0).astype(int)
tb2ea.index = range(101)
tbea = pd.concat([tb2ea, tbea], axis=1)
tbea.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_minpip_vs_lit_ea_stats.tsv", sep='\t')
tbeither = pd.DataFrame(tbeither).fillna(0).astype(int)
tbeither.index = range(101)
tbeither.columns = ["[0,0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.5]", "(0.5,0.9]","(0.9,1]", "1"]
tb2either = pd.DataFrame(tb2either).fillna(0).astype(int)
tb2either.index = range(101)
tbeither = pd.concat([tb2either, tbeither], axis=1)
tbeither.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_minpip_vs_lit_either_stats.tsv", sep='\t')

#effect size comparison:
#annotate beta on aa and ea
aa["beta_aa"] = aa.iloc[:,9]
aa.loc[aa.iloc[:,4]!=aa.iloc[:,5],"beta_aa"] = aa.loc[aa.iloc[:,4]!=aa.iloc[:,5],"beta_aa"]*-1
ea["beta_ea"] = ea.iloc[:,9]
ea.loc[ea.iloc[:,4]!=ea.iloc[:,5],"beta_ea"] = ea.loc[ea.iloc[:,4]!=ea.iloc[:,5],"beta_ea"]*-1
#Filter to those that matter: PIP>0.01 evidence
aasignif = aa[aa.iloc[:,14]>0.01]
easignif = ea[ea.iloc[:,14]>0.01]
pjsignif = pj[pj.p_min_pip>0.01]
#annotate the marginal beta for these pj:
jtsignif = []
for chr in list(range(1,22+1))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(chr)
    df = pd.read_csv(input_file, sep='\t')
    df["variant_id_hg38"] = df.variant_id.str.split("_").str[0]
    df.set_index(["variant_id_hg38","gene_id"], inplace=True)
    dfkeep = df.loc[df.index.intersection(pjsignif.index), :]
    jtsignif.append(dfkeep)
    print("done chr{0}, {1}".format(chr, tm.ctime()))
    print("kept {0} out of {1} lines".format(dfkeep.shape[0], df.shape[0]))
jtsignif = pd.concat(jtsignif)
#take the intersection with aa:
jtsignif.index.names = ['variant_hg38', 'ensg_id']
pjsignif = pjsignif.join(jtsignif.slope, how="left")
aasignif.set_index(["variant_hg38", "ensg_id"], inplace=True)
pjsignif = pjsignif.join(aasignif.loc[:,["14","beta_aa"]], how="left")
pjsignif.rename(columns = {"14":"pip_aa"},inplace=True)
pjsignif = pjsignif.join(easignif.loc[:,["14","beta_ea"]], how="left")
pjsignif.rename(columns = {"14":"pip_ea"},inplace=True)
pjsignif.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eff_size_comparison_vs_ea_aa_df.tsv", sep="\t")

#vep annotation for everything, larger than PIP>0.1 in any category
p_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_finemap_pip0001.txt", sep=' ')
p_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_susie_pip0001.txt", sep='\t')
#actually we only need PIP>=0.1 is fine for vep
p_fm = p_fm[p_fm.prob>=0.1]
p_fm["variant_id_hg38"] = p_fm.rsid.str.split("_").str[0]
p_fm["gene_id"] = p_fm.rsid.str.split("_").str[-1]
p_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
p_sus = p_sus[p_sus.pip>=0.1]
p_sus["variant_id_hg38"] = p_sus.rsid.str.split("_").str[0]
p_sus["gene_id"] = p_sus.rsid.str.split("_").str[-1]
p_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
pj = p_fm.join(p_sus.pip, how="outer")
pj.fillna(0, inplace=True) #anything lower than 0001 is 0
pj.columns = ["rsid", "p_pip_fm", "p_pip_sus"]
pj["p_min_pip"] = np.minimum(pj.p_pip_fm, pj.p_pip_sus)
pj = pj[pj["p_min_pip"]>0.1]
pj.reset_index(inplace=True)
pj["Gene"] = pj.gene_id.str.split("\\.").str[0]
pj.set_index(["variant_id_hg38", "Gene"], inplace=True)
#and zhang et al
ensgids = pd.read_csv("~/Desktop/taskforce_n1102/gencode.v30.genes.parsed.tsv", sep='\t')
ensgids.index =ensgids.gene_name
aa = pd.read_csv("~/Desktop/resources/pqtl_aa_all.tsv.gz", sep="\t")
aa.index = aa.gene_name
aa = aa.join(ensgids.ensg_id, how="left")
aa = aa[aa.iloc[:,14]>0.1]
ea = pd.read_csv("~/Desktop/resources/pqtl_ea_all.tsv.gz", sep="\t")
ea.index = ea.gene_name
ea = ea.join(ensgids.ensg_id, how="left")
ea = ea[ea.iloc[:,14]>0.1]
aa["variant_hg38"] = "chr"+aa.iloc[:,0].astype(str)+":"+aa.iloc[:,1].astype(str)+":"+aa.iloc[:,3]+":"+aa.iloc[:,4]
ea["variant_hg38"] = "chr"+ea.iloc[:,0].astype(str)+":"+ea.iloc[:,1].astype(str)+":"+ea.iloc[:,3]+":"+ea.iloc[:,4]

#vep the unions:
tovep = pd.DataFrame(np.union1d(np.union1d(aa.variant_hg38, ea.variant_hg38), pj.index.get_level_values(0)))
tovep.columns = ["ID"]
tovep["chr"] = tovep.ID.str.split(":").str[0].str.replace("chr","")
tovep["pos"] = tovep.ID.str.split(":").str[1].astype(int)
tovep["ref"] = tovep.ID.str.split(":").str[2]
tovep["alt"] = tovep.ID.str.split(":").str[3]
tovep[["chr","pos", "ID", "ref", "alt"]].sort_values(by=["chr", "pos"]).to_csv("~/Desktop/taskforce_n1102/n1300/basic_stats/pqtl_ea_aa_tovep_sorted.vcf", sep="\t", index=False, header=None)

#and get back the vep result:
vps = pd.read_csv("~/Desktop/taskforce_n1102/n1300/basic_stats/pqtl_jctf_ea_aa_01_vepped.txt", sep='\t')
vps["variant_hg38"] = vps["#Uploaded_variation"]
consc = vps.groupby(["variant_hg38","Gene"]).Consequence.apply(lambda x: ",".join(list(set(x))))
tfb = vps[vps.Gene=="-"].groupby(["variant_hg38"]).Consequence.apply(lambda x: ",".join(list(set(x)))) #per variant consequence, focusing on non-coding ones
vps[vps.Gene!="-"].groupby(["variant_hg38"]).Consequence.apply(lambda x: ",".join(list(set(x)))).value_counts() #per variant consequence, focusing on non-coding ones
#annotate AA
aa["Gene"] = aa.ensg_id.str.split("\\.").str[0]
aa.set_index(["variant_hg38", "Gene"], inplace=True)
aa = aa.join(consc, how="left")
aa.Consequence.fillna("NA", inplace=True)
aa = aa.join(tfb, how="left", rsuffix='_nc', on="variant_hg38")
aa.Consequence_nc.fillna("NA", inplace=True)

#do the "most severe annot" ish-thing as we did before:
aa["causal_class"] = "Others/Unknown"
aa.loc[aa.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
aa.loc[aa.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
aa.loc[aa.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
aa.loc[aa.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
aa.loc[aa.Consequence.str.contains("Missense|missense"), "causal_class"] = "Missense"
aa.loc[aa.Consequence.str.contains('frameshift_variant|stop_gained'), "causal_class"] = "PTV"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = aa[aa.iloc[:,14]>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_updated.tsv", sep="\t")

#EA
ea["Gene"] = ea.ensg_id.str.split("\\.").str[0]
ea.set_index(["variant_hg38", "Gene"], inplace=True)
ea = ea.join(consc, how="left")
ea.Consequence.fillna("NA", inplace=True)
ea = ea.join(tfb, how="left", rsuffix='_nc', on="variant_hg38")
ea.Consequence_nc.fillna("NA", inplace=True)
ea["causal_class"] = "Others/Unknown"
ea.loc[ea.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
ea.loc[ea.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
ea.loc[ea.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
ea.loc[ea.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
ea.loc[ea.Consequence.str.contains("Missense|missense"), "causal_class"] = "Missense"
ea.loc[ea.Consequence.str.contains('frameshift_variant|stop_gained'), "causal_class"] = "PTV"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = ea[ea.iloc[:,14]>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_updated.tsv", sep="\t")

#and our JCTF:
pj.index.names = consc.index.names
pj = pj.join(consc, how="left")
pj.Consequence.fillna("NA", inplace=True)
pj = pj.join(tfb, how="left", rsuffix='_nc', on="variant_hg38")
pj.Consequence_nc.fillna("NA", inplace=True)
pj["causal_class"] = "Others/Unknown"
pj.loc[pj.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
pj.loc[pj.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
pj.loc[pj.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
pj.loc[pj.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
pj.loc[pj.Consequence.str.contains("Missense|missense"), "causal_class"] = "Missense"
pj.loc[pj.Consequence.str.contains('frameshift_variant|stop_gained'), "causal_class"] = "PTV"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = pj[pj.p_min_pip>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_annots_updated.tsv", sep="\t")

#e.g. filter out coding guys to begin with
pjnc = pj.loc[~pj.Consequence.str.contains("synonymous|missense|frameshift_variant|stop_gained")]
pjnc["causal_class"] = "Others/Unknown"
pjnc.loc[pjnc.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
pjnc.loc[pjnc.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
pjnc.loc[pjnc.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
pjnc.loc[pjnc.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = pjnc[pjnc.p_min_pip>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_annots_nc.tsv", sep="\t")

eanc = ea.loc[~ea.Consequence.str.contains("synonymous|missense|frameshift_variant|stop_gained")]
eanc["causal_class"] = "Others/Unknown"
eanc.loc[eanc.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
eanc.loc[eanc.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
eanc.loc[eanc.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
eanc.loc[eanc.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = eanc[eanc.iloc[:,14]>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_nc.tsv", sep="\t")

aanc = aa.loc[~aa.Consequence.str.contains("synonymous|missense|frameshift_variant|stop_gained")]
aanc["causal_class"] = "Others/Unknown"
aanc.loc[aanc.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
aanc.loc[aanc.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
aanc.loc[aanc.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
aanc.loc[aanc.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = aanc[aanc.iloc[:,14]>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_nc.tsv", sep="\t")

#And when filtering missense/PTV and keeping synonymous
pjnc = pj.loc[~pj.Consequence.str.contains("missense|frameshift_variant|stop_gained")]
pjnc["causal_class"] = "Others/Unknown"
pjnc.loc[pjnc.Consequence.str.contains("synonymous"), "causal_class"] = "Synonymous"
pjnc.loc[pjnc.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
pjnc.loc[pjnc.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
pjnc.loc[pjnc.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
pjnc.loc[pjnc.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = pjnc[pjnc.p_min_pip>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_annots_nonseqchange.tsv", sep="\t")

eanc = ea.loc[~ea.Consequence.str.contains("missense|frameshift_variant|stop_gained")]
eanc["causal_class"] = "Others/Unknown"
eanc.loc[eanc.Consequence.str.contains("synonymous"), "causal_class"] = "Synonymous"
eanc.loc[eanc.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
eanc.loc[eanc.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
eanc.loc[eanc.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
eanc.loc[eanc.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = eanc[eanc.iloc[:,14]>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_nonseqchange.tsv", sep="\t")

aanc = aa.loc[~aa.Consequence.str.contains("missense|frameshift_variant|stop_gained")]
aanc["causal_class"] = "Others/Unknown"
aanc.loc[aanc.Consequence.str.contains("synonymous"), "causal_class"] = "Synonymous"
aanc.loc[aanc.Consequence_nc.str.contains("TF_binding_site_variant"), "causal_class"] = "TF binding"
aanc.loc[aanc.Consequence.str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
aanc.loc[aanc.Consequence.str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
aanc.loc[aanc.Consequence.str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = aanc[aanc.iloc[:,14]>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_nonseqchange.tsv", sep="\t")


#Filter coding guys and look at baseline annotations:
pjnc = pj.loc[~pj.Consequence.str.contains("synonymous|missense|frameshift_variant|stop_gained")]
eanc = ea.loc[~ea.Consequence.str.contains("synonymous|missense|frameshift_variant|stop_gained")]
aanc = aa.loc[~aa.Consequence.str.contains("synonymous|missense|frameshift_variant|stop_gained")]
pjnc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_noncoding.tsv.gz", sep="\t", compression='gzip')
eanc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_ea_pip01_noncoding.tsv.gz", sep="\t", compression='gzip')
aanc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_aa_pip01_noncoding.tsv.gz", sep="\t", compression='gzip')


#For each of them, check the baseline annotations intersection:
pjnc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_noncoding.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
eanc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_ea_pip01_noncoding.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
aanc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_aa_pip01_noncoding.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])

import pybedtools
pjnc["chr"] = pjnc.index.get_level_values(0).str.split(":").str[0]
pjnc["start"] = pjnc.index.get_level_values(0).str.split(":").str[1].astype(int) -1
pjnc["end"] = pjnc.start + 1 #assuming SNVs, which is trivial for a missense variant
variants_bed = pjnc.iloc[:,-3:]
variants_bed.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/jctf_pqtl_pip01_noncoding.bed", sep='\t', index=False, header=False)
eanc["chr"] = eanc.index.get_level_values(0).str.split(":").str[0]
eanc["start"] = eanc.index.get_level_values(0).str.split(":").str[1].astype(int) -1
eanc["end"] = eanc.start + 1 #assuming SNVs, which is trivial for a missense variant
variants_bed = eanc.iloc[:,-3:]
variants_bed.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/aric_ea_pqtl_pip01_noncoding.bed", sep='\t', index=False, header=False)
aanc["chr"] = aanc.index.get_level_values(0).str.split(":").str[0]
aanc["start"] = aanc.index.get_level_values(0).str.split(":").str[1].astype(int) -1
aanc["end"] = aanc.start + 1 #assuming SNVs, which is trivial for a missense variant
variants_bed = aanc.iloc[:,-3:]
variants_bed.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/aric_aa_pqtl_pip01_noncoding.bed", sep='\t', index=False, header=False)
#write the bed files intersection:
fns = ['/Users/qingbowang/Desktop/resources/baseline/hg38/DHS_Trynka.bed.hg38.bed',
       '/Users/qingbowang/Desktop/resources/baseline/hg38/DHS_peaks_Trynka.bed.hg38.bed',
       '/Users/qingbowang/Desktop/resources/baseline/hg38/WeakEnhancer_Hoffman.bed.hg38.bed',
       '/Users/qingbowang/Desktop/resources/baseline/hg38/Enhancer_Hoffman.bed.hg38.bed',
       '/Users/qingbowang/Desktop/resources/baseline/hg38/SuperEnhancer_Hnisz.bed.hg38.bed',
       '/Users/qingbowang/Desktop/resources/baseline/hg38/PromoterFlanking_Hoffman.bed.hg38.bed',
       '/Users/qingbowang/Desktop/resources/baseline/hg38/Promoter_UCSC.bed.hg38.bed']
categnames = pd.Series(fns).apply(lambda x: x.replace("/Users/qingbowang/Desktop/resources/baseline/hg38/", "").replace(".bed.hg38.bed",""))
for dt in ["jctf", "aric_aa", "aric_ea"]:
    fn = "/Users/qingbowang/Desktop/taskforce_n1102/tmp/{0}_pqtl_pip01_noncoding.bed".format(dt)
    qbed = pybedtools.BedTool(fn)
    for i in range(len(fns)):
        fn = fns[i]
        categ = categnames[i]
        feat = pybedtools.BedTool(fn)
        it = (qbed + feat) #if intersection exists, the variant survives
        fnout = "/Users/qingbowang/Desktop/taskforce_n1102/tmp/{0}_pqtl_pip01_noncoding_{1}_intsct.bed".format(dt, categ)
        it.saveas(fnout)#write the intersection as bed
#annotate back:
pjnc["variant_pos"] = pjnc.index.get_level_values(0).str.split(":").str[0]+"_"+pjnc.index.get_level_values(0).str.split(":").str[1]
pjnc.set_index("variant_pos", append=True, inplace=True)
for i in range(len(fns)):
    fn = fns[i]
    categ = categnames[i]
    fnintsct = "/Users/qingbowang/Desktop/taskforce_n1102/tmp/jctf_pqtl_pip01_noncoding_{0}_intsct.bed".format(categ)
    try:
        intsct = pd.read_csv(fnintsct, sep="\t", header=None)
        intsct[categ] = True
        intsct.index = intsct.iloc[:,0]+"_"+intsct.iloc[:,2].astype(str) #end is the variant position
        intsct["idx"] = intsct.index
        intsct.drop_duplicates(subset="idx", inplace=True)
        intsct.index.names = ["variant_pos"]
        if i<4:
            print (pjnc.head())
            print (intsct.head())
        pjnc = pjnc.join(intsct[categ], how="left", on="variant_pos")
        pjnc[categ].fillna(False, inplace=True)
    except:
        pjnc[categ] = False #when there's no matching annotation
eanc["variant_pos"] = eanc.index.get_level_values(0).str.split(":").str[0]+"_"+eanc.index.get_level_values(0).str.split(":").str[1]
eanc.set_index("variant_pos", append=True, inplace=True)
for i in range(len(fns)):
    fn = fns[i]
    categ = categnames[i]
    fnintsct = "/Users/qingbowang/Desktop/taskforce_n1102/tmp/aric_ea_pqtl_pip01_noncoding_{0}_intsct.bed".format(categ)
    try:
        intsct = pd.read_csv(fnintsct, sep="\t", header=None)
        intsct[categ] = True
        intsct.index = intsct.iloc[:,0]+"_"+intsct.iloc[:,2].astype(str) #end is the variant position
        intsct["idx"] = intsct.index
        intsct.drop_duplicates(subset="idx", inplace=True)
        intsct.index.names = ["variant_pos"]
        if i<4:
            print (eanc.head())
            print (intsct.head())
        eanc = eanc.join(intsct[categ], how="left", on="variant_pos")
        eanc[categ].fillna(False, inplace=True)
    except:
        eanc[categ] = False #when there's no matching annotation
aanc["variant_pos"] = aanc.index.get_level_values(0).str.split(":").str[0]+"_"+aanc.index.get_level_values(0).str.split(":").str[1]
aanc.set_index("variant_pos", append=True, inplace=True)
for i in range(len(fns)):
    fn = fns[i]
    categ = categnames[i]
    fnintsct = "/Users/qingbowang/Desktop/taskforce_n1102/tmp/aric_aa_pqtl_pip01_noncoding_{0}_intsct.bed".format(categ)
    try:
        intsct = pd.read_csv(fnintsct, sep="\t", header=None)
        intsct[categ] = True
        intsct.index = intsct.iloc[:,0]+"_"+intsct.iloc[:,2].astype(str) #end is the variant position
        intsct["idx"] = intsct.index
        intsct.drop_duplicates(subset="idx", inplace=True)
        intsct.index.names = ["variant_pos"]
        if i<4:
            print (aanc.head())
            print (intsct.head())
        aanc = aanc.join(intsct[categ], how="left", on="variant_pos")
        aanc[categ].fillna(False, inplace=True)
    except:
        aanc[categ] = False #when there's no matching annotation
#save:
pjnc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip')
eanc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_ea_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip')
aanc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_aa_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip')

#get the stats:
pjnc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
eanc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_ea_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
aanc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_aa_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])

pjnc["causal_class"] = "Others/Unknown"
pjnc.loc[pjnc.Consequence.fillna("NA").str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
pjnc.loc[pjnc.Consequence.fillna("NA").str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
pjnc.loc[pjnc.Consequence.fillna("NA").str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
pjnc.loc[pjnc.DHS_Trynka, "causal_class"] = "DHS"
pjnc.loc[pjnc.Enhancer_Hoffman, "causal_class"] = "Enhancer"
pjnc.loc[pjnc.Promoter_UCSC, "causal_class"] = "Promoter"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = pjnc[pjnc.p_min_pip>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pj_pip01_annots_baseline_table.tsv", sep="\t")
print (tb)

eanc["causal_class"] = "Others/Unknown"
eanc.loc[eanc.Consequence.fillna("NA").str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
eanc.loc[eanc.Consequence.fillna("NA").str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
eanc.loc[eanc.Consequence.fillna("NA").str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
eanc.loc[eanc.DHS_Trynka, "causal_class"] = "DHS"
eanc.loc[eanc.Enhancer_Hoffman, "causal_class"] = "Enhancer"
eanc.loc[eanc.Promoter_UCSC, "causal_class"] = "Promoter"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = eanc[eanc.loc[:,"14"]>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_baseline_table.tsv", sep="\t")
print (tb)

aanc["causal_class"] = "Others/Unknown"
aanc.loc[aanc.Consequence.fillna("NA").str.contains("5_prime_UTR_variant"), "causal_class"] = "5'UTR"
aanc.loc[aanc.Consequence.fillna("NA").str.contains("3_prime_UTR_variant"), "causal_class"] = "3'UTR"
aanc.loc[aanc.Consequence.fillna("NA").str.contains("Splice|splice"), "causal_class"] = "Possibly splice"
aanc.loc[aanc.DHS_Trynka, "causal_class"] = "DHS"
aanc.loc[aanc.Enhancer_Hoffman, "causal_class"] = "Enhancer"
aanc.loc[aanc.Promoter_UCSC, "causal_class"] = "Promoter"
tb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = aanc[aanc.loc[:,"14"]>=i].causal_class.value_counts()
    tb.append(v)
tb = pd.concat(tb, axis=1).fillna(0).astype(int)
tb.columns = np.arange(0.1, 1.1, 0.1)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_baseline_table.tsv", sep="\t")
print (tb)

tbj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pj_pip01_annots_baseline_table.tsv", sep="\t", index_col=0)
tbe = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_baseline_table.tsv", sep="\t", index_col=0)
tba = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_baseline_table.tsv", sep="\t", index_col=0)
#EMS overlay:
gpip = pd.read_csv("/Users/qingbowang/Desktop/resources/Whole_Blood_allpairs_gtex_pip.tsv.gz", sep="\t") #takes some time, large file
gpip["gene"] = gpip.gene_id.str.split("\\.").str[0] #takes some time..
gpip.set_index(["variant_id", "gene"], inplace=True)

pjnc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
eanc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_ea_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
aanc = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_aa_pip01_noncoding_baselineannot.tsv.gz", sep="\t", compression='gzip', index_col=[0,1])
pjnc["gene"] = pjnc.index.get_level_values(1)
eanc["gene"] = eanc.index.get_level_values(1)
aanc["gene"] = aanc.index.get_level_values(1)
pjnc["variant_id"] = pjnc.index.get_level_values(0).str.replace(":","_")+"_b38"
eanc["variant_id"] = eanc.index.get_level_values(0).str.replace(":","_")+"_b38"
aanc["variant_id"] = aanc.index.get_level_values(0).str.replace(":","_")+"_b38"
pjnc.set_index(["variant_id", "gene"], inplace=True)
eanc.set_index(["variant_id", "gene"], inplace=True)
aanc.set_index(["variant_id", "gene"], inplace=True)
pjnc = pjnc.join(gpip.ems_bin, how="left")
eanc = eanc.join(gpip.ems_bin, how="left")
aanc = aanc.join(gpip.ems_bin, how="left")
pjnc = pjnc[~pjnc.ems_bin.isna()]
eanc = eanc[~eanc.ems_bin.isna()]
aanc = aanc[~aanc.ems_bin.isna()]
#save:
pjnc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_noncoding_ems.tsv.gz", sep="\t", compression='gzip')
eanc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_ea_pip01_noncoding_ems.tsv.gz", sep="\t", compression='gzip')
aanc.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aric_aa_pip01_noncoding_ems.tsv.gz", sep="\t", compression='gzip')
#and save table:
pjtb = []
eatb = []
aatb = []
for i in np.arange(0.1, 1.1, 0.1):
    v = pjnc[pjnc.p_min_pip >= i].ems_bin.value_counts()
    pjtb.append(v)
    v = eanc[eanc.loc[:, "14"] >= i].ems_bin.value_counts()
    eatb.append(v)
    v = aanc[aanc.loc[:,"14"]>=i].ems_bin.value_counts()
    aatb.append(v)
pjtb = pd.concat(pjtb, axis=1).fillna(0).astype(int).sort_index()
pjtb.columns = np.arange(0.1, 1.1, 0.1)
eatb = pd.concat(eatb, axis=1).fillna(0).astype(int).sort_index()
eatb.columns = np.arange(0.1, 1.1, 0.1)
aatb = pd.concat(aatb, axis=1).fillna(0).astype(int).sort_index()
aatb.columns = np.arange(0.1, 1.1, 0.1)
pjtb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pj_pip01_ems_table.tsv", sep="\t")
eatb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_ems_table.tsv", sep="\t")
aatb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_ems_table.tsv", sep="\t")
print (pjtb/pjtb.sum())
print (eatb/eatb.sum())
print (aatb/aatb.sum())

#vep annotations for pQTLs overall:
p_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_finemap_pip0001.txt", sep=' ')
p_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_susie_pip0001.txt", sep='\t')
p_fm["variant_id_hg38"] = p_fm.rsid.str.split("_").str[0]
p_fm["gene_id"] = p_fm.rsid.str.split("_").str[-1]
p_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
p_sus["variant_id_hg38"] = p_sus.rsid.str.split("_").str[0]
p_sus["gene_id"] = p_sus.rsid.str.split("_").str[-1]
p_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
pj = p_fm.join(p_sus.pip, how="inner") #MIN>0.001
pj.columns = ["rsid", "p_pip_fm", "p_pip_sus"]
pj["p_min_pip"] = np.minimum(pj.p_pip_fm, pj.p_pip_sus)
pj["min_pip_bin"] = pj.p_min_pip.apply(lambda x: bin_pip(x))
pj.reset_index(inplace=True)
pj["Gene"] = pj.gene_id.str.split("\\.").str[0]
pj.set_index(["variant_id_hg38", "Gene"], inplace=True)
vps = []
for chr in list(range(1,23))+["X"]:
    vp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/all_n1419_vepped/n1419_vepped_chr{0}.txt".format(chr), sep='\t')
    vp["variant_id_hg38"] = vp["#Uploaded_variation"].str.split("_").str[0]
    vp.set_index(["variant_id_hg38","Gene"], inplace=True)
    intsct = vp.index.intersection(pj.index)
    vp = vp.loc[intsct,:]
    vps.append(vp)
    if chr==1:
        print (vp)
    print ("done chr{0}, {1}".format(chr, tm.ctime()))
vps = pd.concat(vps)
vps.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_pqtl_pip0001.tsv.gz", sep='\t', compression='gzip')

vps = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_pqtl_pip0001.tsv.gz", sep='\t', index_col=[0,1])
cons = vps.groupby(["variant_id_hg38","Gene"]).Consequence.apply(lambda x: ','.join(x))
pj = pj.join(cons, how="left")
pj.Consequence.fillna("NA", inplace=True)
#get the enrichment, basic annotations:
uniqueannot = pd.Series((",".join(pj.Consequence.unique())).split(",")).unique()
tb = {}
for annot in uniqueannot:
    pj[annot] = pj.Consequence.str.contains(annot)
    t = pj.groupby(["min_pip_bin", annot]).size().unstack().T.fillna(0).astype(int)
    oth = t.sum(axis=1)
    t = t.T
    t["frac"] = t[True]/t.sum(axis=1)
    t["err"] = np.sqrt(t.frac*(1-t.frac)/t.sum(axis=1))
    t.sort_index(inplace=True)
    tb[annot] = t
    print (t)
    t.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pqtl_enr_{0}_updated.tsv".format(annot), sep='\t')
#add two manual combined annot
pj["inframe"] = pj.Consequence.str.contains("inframe_deletion|inframe_insertion|start_lost|stop_lost")
pj["lof"] = pj.Consequence.str.contains("frameshift_variant|stop_gained|splice_acceptor_variant|splice_donor_variant")
for annot in ["inframe", "lof"]:
    t = pj.groupby(["min_pip_bin", annot]).size().unstack().T.fillna(0).astype(int)
    oth = t.sum(axis=1)
    t = t.T
    t["frac"] = t[True]/t.sum(axis=1)
    t["err"] = np.sqrt(t.frac*(1-t.frac)/t.sum(axis=1))
    t.sort_index(inplace=True)
    tb[annot] = t
    print (t)
    t.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pqtl_enr_{0}_updated.tsv".format(annot), sep='\t')

#and then background:
for chr in list(range(1,23))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(chr)
    p = pd.read_csv(input_file, sep="\t")
    p["variant_id_hg38"] = p.variant_id.str.split("_").str[0]
    p["Gene"] = p.gene_id.str.split("\\.").str[0]
    p.set_index(["variant_id_hg38","Gene"], inplace=True)
    vp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/all_n1419_vepped/n1419_vepped_chr{0}.txt".format(chr),sep='\t')
    vp["variant_id_hg38"] = vp["#Uploaded_variation"].str.split("_").str[0]
    vp.set_index(["variant_id_hg38","Gene"], inplace=True)
    cons = vp.groupby(["variant_id_hg38", "Gene"]).Consequence.apply(lambda x: ','.join(x))
    p = p.join(cons, how="left")
    p.Consequence.fillna("NA", inplace=True)
    tb = {}
    for annot in uniqueannot:
        t = p.Consequence.str.contains(annot).value_counts()
        tb[annot] = t
    #my combined annot as well:
    p["inframe"] = p.Consequence.str.contains("inframe_deletion|inframe_insertion|start_lost|stop_lost")
    p["lof"] = p.Consequence.str.contains("frameshift_variant|stop_gained|splice_acceptor_variant|splice_donor_variant")
    for annot in ["inframe","lof"]:
        t = p[annot].value_counts()
        tb[annot] = t
    tb = pd.DataFrame(tb)
    tb.fillna(0).astype(int).to_csv("/Users/qingbowang/Desktop/tmp/vepannot_pqtl_baseline_stats_chr{0}_updated.tsv".format(chr), sep='\t')
    print ("done chr{0}, {1}".format(chr, tm.ctime()))

tb = pd.read_csv("/Users/qingbowang/Desktop/tmp/vepannot_pqtl_baseline_stats_chr{0}_updated.tsv".format(1), sep='\t', index_col=0)
for chr in list(range(2,23))+["X"]:
    tb = tb.add(pd.read_csv("/Users/qingbowang/Desktop/tmp/vepannot_pqtl_baseline_stats_chr{0}_updated.tsv".format(chr), sep='\t', index_col=0), fill_value=0)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_pqtl_baseline_stats_updated.tsv".format(chr), sep='\t')
annot_to_use = ["intron_variant", "upstream_gene_variant", "downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant",
              "synonymous_variant", "missense_variant","inframe", "lof"]#do we need this at all? Just outputting everything should be fine

#Then deeper investigation focusing on missense variants:
pj = pj[pj.Consequence.str.contains("missense")]
#save once
pj.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_pip0001_missense_updated.tsv.gz", sep='\t', compression="gzip")
#and annotate various biochemical features:
pj["chr"] = pj.index.get_level_values(0).str.split(":").str[0]
pj["start"] = pj.index.get_level_values(0).str.split(":").str[1].astype(int) -1
pj["end"] = pj.start + 1 #assuming SNVs, which is trivial for a missense variant
variants_bed = pj.iloc[:,-3:]
variants_bed.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_pip0001_aminoacids_updated.bed", sep='\t', index=False, header=False)#these correspond to missense vars of PIP>0.001
#and get the intersection for each bed files:
import pybedtools
fn = "/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_pip0001_aminoacids_updated.bed"
aabed = pybedtools.BedTool(fn)
import glob
fns = glob.glob('/Users/qingbowang/Desktop/resources/uniprot_annots/uniprot_*_clean.bed')
categnames = pd.Series(fns).apply(lambda x: x.replace("/Users/qingbowang/Desktop/resources/uniprot_annots/uniprot_", "").replace("_clean.bed",""))
for i in range(len(fns)):
    fn = fns[i]
    categ = categnames[i]
    feat = pybedtools.BedTool(fn)
    it = (aabed + feat) #if intersection exists, the variant survives
    fnout = "/Users/qingbowang/Desktop/resources/uniprot_annots/out/uniprot_{0}_pip0001_intsct_updated.bed".format(categ)
    it.saveas(fnout)#write the intersection as bed
#annotate back the original file:
pj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_pip0001_missense_updated.tsv.gz", sep='\t')
pj.index = pj.variant_id_hg38.str.split(":").str[0]+"_"+ pj.variant_id_hg38.str.split(":").str[1]
for i in range(len(fns)):
    fn = fns[i]
    categ = categnames[i]
    fnintsct = "/Users/qingbowang/Desktop/resources/uniprot_annots/out/uniprot_{0}_pip0001_intsct_updated.bed".format(categ)
    try:
        intsct = pd.read_csv(fnintsct, sep="\t", header=None)
        intsct[categ] = True
        intsct.index = intsct.iloc[:,0]+"_"+intsct.iloc[:,2].astype(str) #end is the variant position
        intsct["idx"] = intsct.index
        intsct.drop_duplicates(subset="idx", inplace=True)
        if i<4:
            print (pj.head())
            print (intsct.head())
        pj = pj.join(intsct[categ], how="left")
        pj[categ].fillna(False, inplace=True)
    except:
        pj[categ] = False #when there's no matching annotation
pj.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/vepannot_pip0001_aminoacids_uniprotannot_updated.tsv.gz", sep='\t', index=False, compression="gzip")

#Also background
tbs = []
for chr in list(range(1,23))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(str(chr).replace("head","").replace("tail",""))
    p = pd.read_csv(input_file, sep="\t")
    p["variant_id_hg38"] = p.variant_id.str.split("_").str[0]
    p["Gene"] = p.gene_id.str.split("\\.").str[0]
    p.set_index(["variant_id_hg38","Gene"], inplace=True)
    vp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/all_n1419_vepped/n1419_vepped_chr{0}.txt".format(chr),sep='\t')
    vp["variant_id_hg38"] = vp["#Uploaded_variation"].str.split("_").str[0]
    vp.set_index(["variant_id_hg38","Gene"], inplace=True)
    cons = vp.groupby(["variant_id_hg38", "Gene"]).Consequence.apply(lambda x: ','.join(x))
    p = p.join(cons[cons.str.contains("missense")], how="left")
    tbs.append(p)
    print ("done chr{0}, {1}".format(chr, tm.ctime()))
tbs = pd.concat(tbs)
tbs = tbs[tbs.Consequence.fillna("NA").str.contains("missense")]
tbs.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pqtl_all_aminoacids_updated.tsv", sep='\t')

#and make the bed:
tbs = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pqtl_all_aminoacids_updated.tsv", sep='\t', index_col=[0,1])
tbs["chr"] = tbs.index.get_level_values(0).str.split(":").str[0]
tbs["start"] = tbs.index.get_level_values(0).str.split(":").str[1].astype(int) -1
tbs["end"] = tbs.start + 1 #assuming SNVs, which is trivial for a missense variant
variants_bed = tbs.iloc[:,-3:]
variants_bed.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_all_aminoacids_updated.bed", sep='\t', index=False, header=False)#these correspond to missense vars of PIP>0.001
#and get the intersection for each bed files:
fn = "/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_all_aminoacids_updated.bed"
aabed = pybedtools.BedTool(fn)
import glob
fns = glob.glob('/Users/qingbowang/Desktop/resources/uniprot_annots/uniprot_*_clean.bed')
categnames = pd.Series(fns).apply(lambda x: x.replace("/Users/qingbowang/Desktop/resources/uniprot_annots/uniprot_", "").replace("_clean.bed",""))
for i in range(len(fns)):
    fn = fns[i]
    categ = categnames[i]
    feat = pybedtools.BedTool(fn)
    it = (aabed + feat) #if intersection exists, the variant survives
    fnout = "/Users/qingbowang/Desktop/resources/uniprot_annots/out/uniprot_{0}_background_intsct_updated.bed".format(categ)
    it.saveas(fnout)#write the intersection as bed
#and put it back to the original data
tbs = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/pqtl_all_aminoacids_updated.tsv", sep='\t')
tbs.index = tbs.variant_id.str.split(":").str[0]+"_"+ tbs.variant_id.str.split(":").str[1]
for i in range(len(fns)):
    fn = fns[i]
    categ = categnames[i]
    fnintsct = "/Users/qingbowang/Desktop/resources/uniprot_annots/out/uniprot_{0}_background_intsct_updated.bed".format(categ)
    try:
        intsct = pd.read_csv(fnintsct, sep="\t", header=None)
        intsct[categ] = True
        intsct.index = intsct.iloc[:,0]+"_"+intsct.iloc[:,2].astype(str) #end is the variant position
        intsct["idx"] = intsct.index
        intsct.drop_duplicates(subset="idx", inplace=True)
        if i<4:
            print (tbs.head())
            print (intsct.head())
        tbs = tbs.join(intsct[categ], how="left")
        tbs[categ].fillna(False, inplace=True)
    except:
        tbs[categ] = False #when there's no matching annotation
tbs.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/vepannot_background_aminoacids_uniprotannot_updated.tsv.gz", sep='\t', index=False, compression="gzip")

#get the summary table:
tbs = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/vepannot_background_aminoacids_uniprotannot_updated.tsv.gz", sep='\t')
interest = ["chains", "structure_beta", "structure_helix", "disulfbonds", "extracellular", "transmembrane", 'cytoplasmic']
tbout = {}
for categ in interest:
    tb = tbs[categ].value_counts()
    tb["frac"] = tb[True] / tb.sum()
    tb["err"] = np.sqrt(tb.frac*(1-tb.frac)/(tb[True]+tb[False]))
    tbout[categ] = tb
    print (tb)
tbout = pd.DataFrame(tbout)
tbout.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_background_uniprot_enr_updated.tsv", sep='\t')

interest = ["chains", "structure_beta", "structure_helix", "disulfbonds", "extracellular", "transmembrane", 'cytoplasmic']
tbout = {}
for categ in interest:
    tb = pj.groupby(["min_pip_bin", categ]).size().unstack().fillna(0).astype(int)
    if tb.shape[1]==1:
        tb[True] = 0
    tb["frac"] = tb[True] / tb.sum(axis=1)
    tb["err"] = np.sqrt(tb.frac*(1-tb.frac)/(tb[True]+tb[False]))
    tbout[categ] = tb
    print (tb)
    tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_uniprot_{0}_enr_updated.tsv".format(categ), sep='\t')

#For non-missense guys, tissue origin analysis
#get the GTEx PIP: (not restricting to PIP>0.1, but let's do all PIP>0.0001) - takes some time since this is the full data.
gtex = pd.read_csv("/Users/qingbowang/Desktop/enformer_prep/susie_results-GTEx_49tissues_release1.tsv.bgz", sep="\t", compression="gzip")
gtex_sus = gtex[gtex["method"]!="FINEMAP"][["variant_hg38","gene", "tissue", "pip"]]
gtex_fm = gtex[gtex["method"]=="FINEMAP"][["variant_hg38","gene", "tissue", "pip"]]
gtex_sus.set_index(["variant_hg38","gene", "tissue"], inplace=True)
gtex_fm.set_index(["variant_hg38","gene", "tissue"], inplace=True)
gtex_sus.columns = ["gtex_pip_susie"]
gtex_fm.columns = ["gtex_pip_fm"]
gtex_sus = gtex_sus.join(gtex_fm.gtex_pip_fm, how="outer").fillna(0)
gtex_sus["gtex_min_pip"] = np.minimum(gtex_sus.gtex_pip_susie, gtex_sus.gtex_pip_fm)
gtex_sus = gtex_sus[gtex_sus.gtex_min_pip>0.0001]
gtex_sus.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/GTEx_49tissues_pip0001.tsv.gz", sep="\t", compression="gzip")

#and overlay the pQTL PIP info:
p_fm = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_finemap_pip0001.txt", sep=' ')
p_sus = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/n1384_pqtl_susie_pip0001.txt", sep='\t')
p_fm["variant_id_hg38"] = p_fm.rsid.str.split("_").str[0]
p_fm["gene_id"] = p_fm.rsid.str.split("_").str[-1]
p_fm.set_index(["variant_id_hg38","gene_id"], inplace=True)
p_sus["variant_id_hg38"] = p_sus.rsid.str.split("_").str[0]
p_sus["gene_id"] = p_sus.rsid.str.split("_").str[-1]
p_sus.set_index(["variant_id_hg38","gene_id"], inplace=True)
pj = p_fm.join(p_sus.pip, how="inner") #MIN>0.001
pj.columns = ["rsid", "p_pip_fm", "p_pip_sus"]
pj["p_min_pip"] = np.minimum(pj.p_pip_fm, pj.p_pip_sus)
pj["min_pip_bin"] = pj.p_min_pip.apply(lambda x: bin_pip(x))
#and exclude missense and lof guys:
vps = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_pqtl_pip0001.tsv.gz", sep='\t', index_col=[0,1])
cons = vps.groupby(["variant_id_hg38","Gene"]).Consequence.apply(lambda x: ','.join(x))
mis = cons[cons.str.contains("missense|frameshift_variant|stop_gained|inframe_deletion|inframe_insertion|start_lost|stop_lost")]
pj["Gene"] = pj.index.get_level_values(1).str.split("\\.").str[0]
pj.set_index("Gene", append=True, inplace=True)
pj = pj.join(mis, how="left", on=["variant_id_hg38", "Gene"]) #adding "Consequence" column only when it contains what we want to exclude

gtex = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/GTEx_49tissues_pip0001.tsv.gz", sep="\t", compression="gzip")
#set index common:
pj["v"] = pj.index.get_level_values(0)
pj["g"] = pj.index.get_level_values(1).str.split("\\.").str[0]
pj.set_index(["v","g"], inplace=True)
gtex["v"] = gtex.variant_hg38.str.replace("_b38","").str.replace("_",":")
gtex["g"] = gtex.gene.str.split("\\.").str[0]
gtex.set_index(["v","g"], inplace=True)
gtex = gtex.join(pj[["p_min_pip","Consequence"]], how="left")
gtex["clpp"] = gtex.p_min_pip*gtex.gtex_min_pip
def bin_clpp(x):
    if x>=0.9: return (0.9)
    elif x>=0.1: return (0.1)
    elif x>=0.01: return (0.01)
    else: return (0)
gtex["clpp_bin"] = gtex.clpp.apply(lambda x: bin_clpp(x))
tb = gtex.groupby(["tissue", "clpp_bin"]).size().unstack().fillna(0).astype(int).iloc[:,1:] #No need for 0 bin
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_table.tsv", sep="\t")

#take the sum CLPP per tissue:
#gtex = gtex[gtex.gtex_min_pip>0.001]#Either <0.001 in one tissue -> 0 no, 0.0001 filter is fine
gtex["clpp"] = gtex["clpp"].fillna(0)
tb_clpp = gtex.groupby("tissue").clpp.sum()
tb_clpp.sort_values()
tb_clpp.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_sum_pertissue.tsv", sep="\t")
n_clpp = gtex.groupby("tissue").size()
n_clpp.sort_values()
n_clpp.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_n_pertissue.tsv", sep="\t")
n_09 = gtex[gtex.gtex_min_pip>0.9].groupby("tissue").size()
n_09.sort_values()
n_09.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_n_pip09_pertissue.tsv", sep="\t")

#Gtex enr:
def test_a_tissue(gtex, iter=1000, tissue="Whole_Blood", filter_missense=True):
    if filter_missense:
        gtex = gtex[gtex.Consequence.isna()]
    dfcase = gtex[gtex.tissue == tissue][["p_min_pip","gtex_min_pip"]]
    N = dfcase.shape[0]
    dfcont = gtex[gtex.tissue != tissue][["p_min_pip","gtex_min_pip"]]
    dfcont["clpp"] = dfcont.p_min_pip*dfcont.gtex_min_pip
    numer = sum(dfcase.p_min_pip*dfcase.gtex_min_pip)
    denom = np.random.choice(dfcont.clpp, size=(iter, N),replace=True).sum(axis=1)
    return ([numer, denom])
vs = {}
i = 0
for tissue in gtex.tissue.unique():
    vs[tissue] = test_a_tissue(gtex, 2000, tissue)
    print ("done {0}, {1}, {2}".format(tissue, i, tm.ctime()))
    i += 1
enr = []
errlow = []
errhigh = []
pval = []
for tissue in gtex.tissue.unique():
    t = vs[tissue]
    enrs = (t[0]/np.array(t[1]))
    enrs.sort()
    enr.append(np.median(enrs))
    errlow.append(enrs[50])
    errhigh.append(enrs[-50])
    pval.append(min(sum(t[0]>np.array(t[1]))/2000, 1 - sum(t[0]>np.array(t[1]))/2000)*2)
tb = pd.DataFrame({"enr":enr, "upper":errhigh, "lower":errlow, "pval":pval})
tb.index = gtex.tissue.unique()
tb.sort_values(by="enr", ascending=False)
tb.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_enrichment_scores.tsv", sep="\t")



### Plotting part
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib.colors import LogNorm
import seaborn as sns

#First, the p-value distribution
plt.rcParams.update({'font.size': 16}) #Updated 20230626
minp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/basic_stats/pqtl_minp_pergene.txt", sep="\t", index_col=0, squeeze=True)
cum = pd.Series(np.arange(201)/10).apply(lambda x: sum(minp<10**-x))
N = int(minp.shape[0])
plt.figure(figsize=(6,4.5))
plt.plot(np.arange(len(cum)), cum/len(minp), color="black")
#plt.xticks(np.arange(len(cum))[::10], (np.arange(201)/10)[::-1][::10].astype(int))
plt.xticks(np.arange(len(cum))[::10], (np.arange(201)/10)[::10].astype(int))
plt.ylim([-0.05,1.05])
plt.axhline(0, linestyle="dotted", zorder=-1, color="black", linewidth=0.5)
plt.axhline(1000/len(minp), linestyle="dotted", zorder=-1, color="black", linewidth=0.5)
plt.axhline(2000/len(minp), linestyle="dotted", zorder=-1, color="black", linewidth=0.5)
plt.axhline(1, linestyle="dotted", zorder=-1, color="black", linewidth=0.5)
plt.text(x=130, y=1000/N, s="n=1000", va="top")
plt.text(x=130, y=2000/N, s="n=2000", va="top")
plt.text(x=130, y=1, s="n={0}\n(All tested genes)".format(N), va="top")
plt.xticks(fontsize=13)
plt.ylabel("Fraction of genes", fontsize=17)
plt.xlabel("pQTL -log$_{10}$(p)$\leq$x", fontsize=17)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pergene_minp.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pergene_minp.pdf', bbox_inches='tight', dpi=500)
plt.clf()
plt.rcParams.update({'font.size': 12}) #Putting back 20230626

#SuSiE vs FM:
from matplotlib.colors import LogNorm
import seaborn as sns
tk = ["[0,0.001)","[0.001,0.01)","[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1)","1"]
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary/pqtl_sus_vs_fm_bins.tsv", sep='\t', index_col=0)
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
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_fm_vs_susie.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_fm_vs_susie.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#vs ARIC study:

#1. Our min PIP vs their PIP max(ea, aa) as the main figure
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_minpip_vs_lit_either_stats.tsv", sep='\t', index_col=0)
tb = tb.iloc[:,3:] #remove missing guys
c7 = ["lightskyblue", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
plt.rcParams.update({'font.size': 14}) #Updated 20230626
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(9, 4.5), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("Count", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[0].set_yticks([10**2, 10**4], fontsize=13)
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c7[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c7[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in ARIC study\n(Zhang et al.):", fontsize=14, title_fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
ax[0].set_ylim([10**0.5,10**5.5]) #manual
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aric_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aric_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()


#1.5 same including NAs as a sub figure
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_minpip_vs_lit_either_stats.tsv", sep='\t', index_col=0)
c10 = ["#bbbbbbff", "#999999ff", "#555555ff", "lightskyblue", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c10[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c10[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in ARIC study\n(Zhang et al.):", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
#ax[0].set_ylim([10**0.5,10**4.5]) #manual
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aric_full_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aric_full_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#1.75 same for ea and aa as sub figures
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_minpip_vs_lit_aa_stats.tsv", sep='\t', index_col=0)
tb = tb.iloc[:,3:] #remove missing guys
c7 = ["lightskyblue", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c7[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c7[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in ARIC study\n(African Ancestry):", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
#ax[0].set_ylim([10**0.5,10**4.5]) #manual
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aric_aa_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aric_aa_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_minpip_vs_lit_ea_stats.tsv", sep='\t', index_col=0)
tb = tb.iloc[:,3:] #remove missing guys
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c7[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c7[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in ARIC study\n(European Ancestry):", fontsize=14)
ax[1].set_xlabel("Significance threshold (PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
#ax[0].set_ylim([10**0.5,10**4.5]) #manual
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aric_ea_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aric_ea_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()


#2. their PIP vs our PIP, as a sub figure
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/lit_aa_vs_our_pqtl_minpip_stats.tsv", sep='\t', index_col=0)
tb = tb.iloc[:,:-3]
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c7[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c7[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in our study:", fontsize=14)
ax[1].set_xlabel("PIP threshold in Zhang et al., African Ancestry", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/aric_aa_vs_us_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/aric_aa_vs_us_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()
#EA
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/lit_ea_vs_our_pqtl_minpip_stats.tsv", sep='\t', index_col=0)
tb = tb.iloc[:,:-3]
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c7[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c7[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in our study:", fontsize=14)
ax[1].set_xlabel("PIP threshold in Zhang et al., European Ancestry", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/aric_ea_vs_us_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/aric_ea_vs_us_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#When excluding NAs, our SuSiE PIP instead of min PIP as sub figures
tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/our_pqtl_suspip_vs_lit_either_stats.tsv", sep='\t', index_col=0)
tb = tb.iloc[:,3:] #remove missing guys
tbnorm = tb.apply(lambda x: x/sum(x), axis=1)
#plot:
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb.index, tb.sum(axis=1), log=True, color="black")
ax[0].set_ylabel("N(variant-gene)", fontsize=15)
yM = tb.sum(axis=1)[0]
ax[0].set_xlim([-1,101])
ax[1].bar(tbnorm.index, tbnorm.iloc[:,0], label=tbnorm.columns[0], color=c7[0])
for i in range(1,tbnorm.shape[1]):
    ax[1].bar(tbnorm.index, tbnorm.iloc[:,i], bottom=tbnorm.iloc[:, :i].sum(axis=1), label=tbnorm.columns[i], color=c7[i])
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="PIP in ARIC study\n(Zhang et al.):", fontsize=14)
ax[1].set_xlabel("Significance threshold (SuSiE PIP)", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_xlim([-1,101])
ax[1].set_ylim([0,1])
#ax[0].set_ylim([10**0.5,10**4.5]) #manual
plt.subplots_adjust(hspace = .05)
ax[1].set_xticklabels(list(np.array(ax[1].get_xticks().tolist())/100))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_susiepip_vs_aric_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_susiepip_vs_aric_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#effect size comparison:
pjsignif = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/eff_size_comparison_vs_ea_aa_df.tsv", sep="\t", index_col=[0,1])
from scipy import stats
#annotate min as min between two studies of comparison
pjsignif["ea_min_pip"] = np.minimum(pjsignif.p_min_pip, pjsignif.pip_ea.fillna(0))
pjsignif["aa_min_pip"] = np.minimum(pjsignif.p_min_pip, pjsignif.pip_aa.fillna(0))
#and plot
plt.figure(figsize=(4,4))
d = pjsignif[pjsignif.aa_min_pip>=0.01]
plt.scatter(d.slope, d.beta_aa, color="tab:blue", label="0.01<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
d = pjsignif[pjsignif.aa_min_pip>=0.1]
plt.scatter(d.slope, d.beta_aa, color="tab:green", label="0.1<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
d = pjsignif[pjsignif.aa_min_pip>=0.9]
plt.scatter(d.slope, d.beta_aa, color="tab:orange", label="0.9<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
plt.plot([-2,2], [-2,2], linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.axvline(x=0, linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.axhline(y=0, linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.title("Comparison with Zhang et al. AFR\n(Pearson r={0})".format(np.round(stats.pearsonr(d.slope, d.beta_aa)[0], 4)), fontsize=14)
plt.xlabel("Effect size in JCTF (our study)", fontsize=14)
plt.ylabel("Effect size in AFR", fontsize=14)
plt.legend(framealpha=0.5, fontsize=11)
plt.gca().set_aspect('equal')
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/us_vs_aa_eff_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/us_vs_aa_eff_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

plt.figure(figsize=(4,4))
d = pjsignif[pjsignif.ea_min_pip>=0.01]
plt.scatter(d.slope, d.beta_ea, color="tab:blue", label="0.01<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
d = pjsignif[pjsignif.ea_min_pip>=0.1]
plt.scatter(d.slope, d.beta_ea, color="tab:green", label="0.1<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
d = pjsignif[pjsignif.ea_min_pip>=0.9]
plt.scatter(d.slope, d.beta_ea, color="tab:orange", label="0.9<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
plt.plot([-2,2], [-2,2], linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.axvline(x=0, linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.axhline(y=0, linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.title("Comparison with Zhang et al. EUR\n(Pearson r={0})".format(np.round(stats.pearsonr(d.slope, d.beta_ea)[0], 4)), fontsize=14)
plt.xlabel("Effect size in JCTF (our study)", fontsize=14)
plt.ylabel("Effect size in EUR", fontsize=14)
plt.legend(framealpha=0.5, fontsize=11)
plt.gca().set_aspect('equal')
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/us_vs_ea_eff_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/us_vs_ea_eff_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()


#Them vs them?
plt.figure(figsize=(4,4))
d = pjsignif[(pjsignif.ea_min_pip>=0.01)&(pjsignif.aa_min_pip>=0.01)]
plt.scatter(d.beta_ea, d.beta_aa, color="tab:blue", label="0.01<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
d = pjsignif[(pjsignif.ea_min_pip>=0.1)&(pjsignif.aa_min_pip>=0.1)]
plt.scatter(d.beta_ea, d.beta_aa, color="tab:green", label="0.1<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
d = pjsignif[(pjsignif.ea_min_pip>=0.9)&(pjsignif.aa_min_pip>=0.9)]
plt.scatter(d.beta_ea, d.beta_aa, color="tab:orange", label="0.9<PIP\n(n={0})".format(d.shape[0]), linewidth=0.2, edgecolors="#ffffffff")
plt.plot([-2,2], [-2,2], linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.axvline(x=0, linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.axhline(y=0, linestyle="--", color="tab:gray", linewidth=0.5, zorder=-1)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.title("Zhang et al. EUR vs AFR\n(Pearson r={0})".format(np.round(stats.pearsonr(d.beta_ea, d.beta_aa)[0], 4)), fontsize=14)
plt.xlabel("Effect size in EUR", fontsize=14)
plt.ylabel("Effect size in AFR", fontsize=14)
plt.legend()
plt.gca().set_aspect('equal')
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/ea_vs_aa_eff_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/ea_vs_aa_eff_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#functional annotations comparison:
tb1 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_annots_updated.tsv", sep="\t", index_col=0)
tb2 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_updated.tsv", sep="\t", index_col=0)
tb3 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_updated.tsv", sep="\t", index_col=0)
tb1.columns = np.arange(0.1, 1.1, 0.1)
tb2.columns = np.arange(0.1, 1.1, 0.1)
tb3.columns = np.arange(0.1, 1.1, 0.1) 
tb1 = tb1.iloc[[6,1,5,3,4,2,0],:]
tb2 = tb2.iloc[[6,1,5,3,4,2,0],:]
tb3 = tb3.iloc[[6,1,5,3,4,2,0],:]
f1 = tb1/tb1.sum()
f2 = tb2/tb2.sum()
f3 = tb3/tb3.sum()
print (f1)
print (f2)
print (f3)
#Plot in a single plot:
c7 = ["tab:red", "tab:orange", "yellow", "tab:green", "yellowgreen", "tab:blue", "tab:gray"]
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(10, 5), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb1.columns-0.03, tb1.sum(axis=0), log=True, width=0.03, color="tab:pink", label="JCTF (our study)")
ax[0].bar(tb3.columns, tb3.sum(axis=0), log=True, width=0.03, color="tab:olive", label="Zhang et al., AFR")
ax[0].bar(tb2.columns+0.03, tb2.sum(axis=0), log=True, width=0.03, color="tab:cyan", label="Zhang et al., EUR")
ax[0].set_ylabel("Count", fontsize=14)
ax[0].set_yticks([10**2, 10**3, 10**4])
ax[0].set_xlim([0.04,1.06])
ax[0].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Study:", fontsize=14, title_fontsize=14)
ax[1].bar(f1.columns-0.03, f1.iloc[0,:], label=f1.index[0], color=c7[0], width=0.03)
for i in range(1,f1.shape[0]):
    ax[1].bar(f1.columns-0.03, f1.iloc[i,:], bottom=f1.iloc[:i,:].sum(axis=0), label=f1.index[i], color=c7[i], width=0.03)
ax[1].bar(f3.columns, f3.iloc[0, :], color=c7[0], width=0.03)
for i in range(1, f3.shape[0]):
    ax[1].bar(f3.columns, f3.iloc[i, :], bottom=f3.iloc[:i, :].sum(axis=0), color=c7[i], width=0.03)
ax[1].bar(f2.columns+0.03, f2.iloc[0,:], color=c7[0], width=0.03)
for i in range(1,f2.shape[0]):
    ax[1].bar(f2.columns+0.03, f2.iloc[i,:], bottom=f2.iloc[:i,:].sum(axis=0), color=c7[i], width=0.03)
#for ax[1], also the edge:
ax[1].bar(f1.columns-0.03, [1]*f1.shape[1], fill = False, edgecolor = "tab:pink", linewidth=2, width=0.028)
ax[1].bar(f3.columns+0.03, [1]*f3.shape[1], fill = False, edgecolor = "tab:cyan", linewidth=2, width=0.028)
ax[1].bar(f2.columns, [1]*f2.shape[1], fill = False, edgecolor = "tab:olive", linewidth=2, width=0.028)
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Variant class:", fontsize=14, title_fontsize=14)
ax[1].set_xlabel("x$\leq$PIP", fontsize=16)
ax[1].set_ylabel("Fraction", fontsize=16)
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticks(np.round(tb1.columns,2))
ax[1].set_xticklabels(np.round(tb1.columns,2))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_vep_vs_aa_and_ea_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_vep_vs_aa_and_ea_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#focusing on non-coding
tb1 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_annots_nc.tsv", sep="\t", index_col=0)
tb2 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_nc.tsv", sep="\t", index_col=0)
tb3 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_nc.tsv", sep="\t", index_col=0)
tb1.columns = np.arange(0.1, 1.1, 0.1)
tb2.columns = np.arange(0.1, 1.1, 0.1)
tb3.columns = np.arange(0.1, 1.1, 0.1) 
tb1 = tb1.iloc[[4,2,3,1,0],:]
tb2 = tb2.iloc[[4,2,3,1,0],:]
tb3 = tb3.iloc[[4,2,3,1,0],:]
f1 = tb1/tb1.sum()
f2 = tb2/tb2.sum()
f3 = tb3/tb3.sum()
print (f1)
print (f2)
print (f3)
#Plot in a single plot:
c5 = ["yellow", "tab:green", "yellowgreen", "tab:blue", "tab:gray"]
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb1.columns-0.03, tb1.sum(axis=0), log=True, width=0.03, color="tab:pink", label="JCTF (our study)")
ax[0].bar(tb3.columns, tb3.sum(axis=0), log=True, width=0.03, color="tab:olive", label="Zhang et al., AFR")
ax[0].bar(tb2.columns+0.03, tb2.sum(axis=0), log=True, width=0.03, color="tab:cyan", label="Zhang et al., EUR")
ax[0].set_ylabel("N(variant-gene)", fontsize=11)
ax[0].set_xlim([0.04,1.06])
ax[0].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Study:", fontsize=14)
ax[1].bar(f1.columns-0.03, f1.iloc[0,:], label=f1.index[0], color=c5[0], width=0.03)
for i in range(1,f1.shape[0]):
    ax[1].bar(f1.columns-0.03, f1.iloc[i,:], bottom=f1.iloc[:i,:].sum(axis=0), label=f1.index[i], color=c5[i], width=0.03)
ax[1].bar(f3.columns, f3.iloc[0, :], color=c5[0], width=0.03)
for i in range(1, f3.shape[0]):
    ax[1].bar(f3.columns, f3.iloc[i, :], bottom=f3.iloc[:i, :].sum(axis=0), color=c5[i], width=0.03)
ax[1].bar(f2.columns+0.03, f2.iloc[0,:], color=c5[0], width=0.03)
for i in range(1,f2.shape[0]):
    ax[1].bar(f2.columns+0.03, f2.iloc[i,:], bottom=f2.iloc[:i,:].sum(axis=0), color=c5[i], width=0.03)
#for ax[1], also the edge:
ax[1].bar(f1.columns-0.03, [1]*f1.shape[1], fill = False, edgecolor = "tab:pink", linewidth=2, width=0.028)
ax[1].bar(f3.columns+0.03, [1]*f3.shape[1], fill = False, edgecolor = "tab:cyan", linewidth=2, width=0.028)
ax[1].bar(f2.columns, [1]*f2.shape[1], fill = False, edgecolor = "tab:olive", linewidth=2, width=0.028)
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Variant class:", fontsize=14)
ax[1].set_xlabel("x$\leq$PIP", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticks(np.round(tb1.columns,2))
ax[1].set_xticklabels(np.round(tb1.columns,2))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_vep_vs_aa_and_ea_noncoding.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aa_and_ea_noncoding.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#focusing on non-coding and synonymous:
tb1 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/jctf_pip01_annots_nonseqchange.tsv", sep="\t", index_col=0)
tb2 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_nonseqchange.tsv", sep="\t", index_col=0)
tb3 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_nonseqchange.tsv", sep="\t", index_col=0)
tb1.columns = np.arange(0.1, 1.1, 0.1)
tb2.columns = np.arange(0.1, 1.1, 0.1)
tb3.columns = np.arange(0.1, 1.1, 0.1)
tb1 = tb1.iloc[[4,2,3,1,5,0],:]
tb2 = tb2.iloc[[4,2,3,1,5,0],:]
tb3 = tb3.iloc[[4,2,3,1,5,0],:]
f1 = tb1/tb1.sum()
f2 = tb2/tb2.sum()
f3 = tb3/tb3.sum()
print (f1)
print (f2)
print (f3)
#Plot in a single plot:
c6 = ["yellow", "tab:green", "yellowgreen", "tab:blue", "lightskyblue", "tab:gray"]
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb1.columns-0.03, tb1.sum(axis=0), log=True, width=0.03, color="tab:pink", label="JCTF (our study)")
ax[0].bar(tb3.columns, tb3.sum(axis=0), log=True, width=0.03, color="tab:olive", label="Zhang et al., AFR")
ax[0].bar(tb2.columns+0.03, tb2.sum(axis=0), log=True, width=0.03, color="tab:cyan", label="Zhang et al., EUR")
ax[0].set_ylabel("N(variant-gene)", fontsize=11)
ax[0].set_xlim([0.04,1.06])
ax[0].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Study:", fontsize=14)
ax[1].bar(f1.columns-0.03, f1.iloc[0,:], label=f1.index[0], color=c6[0], width=0.03)
for i in range(1,f1.shape[0]):
    ax[1].bar(f1.columns-0.03, f1.iloc[i,:], bottom=f1.iloc[:i,:].sum(axis=0), label=f1.index[i], color=c6[i], width=0.03)
ax[1].bar(f3.columns, f3.iloc[0, :], color=c6[0], width=0.03)
for i in range(1, f3.shape[0]):
    ax[1].bar(f3.columns, f3.iloc[i, :], bottom=f3.iloc[:i, :].sum(axis=0), color=c6[i], width=0.03)
ax[1].bar(f2.columns+0.03, f2.iloc[0,:], color=c6[0], width=0.03)
for i in range(1,f2.shape[0]):
    ax[1].bar(f2.columns+0.03, f2.iloc[i,:], bottom=f2.iloc[:i,:].sum(axis=0), color=c6[i], width=0.03)
#for ax[1], also the edge:
ax[1].bar(f1.columns-0.03, [1]*f1.shape[1], fill = False, edgecolor = "tab:pink", linewidth=2, width=0.028)
ax[1].bar(f3.columns+0.03, [1]*f3.shape[1], fill = False, edgecolor = "tab:cyan", linewidth=2, width=0.028)
ax[1].bar(f2.columns, [1]*f2.shape[1], fill = False, edgecolor = "tab:olive", linewidth=2, width=0.028)
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Variant class:", fontsize=14)
ax[1].set_xlabel("x$\leq$PIP", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticks(np.round(tb1.columns,2))
ax[1].set_xticklabels(np.round(tb1.columns,2))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_vep_vs_aa_and_ea_syn.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aa_and_ea_syn.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Focusing on non-coding annot
tb1 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pj_pip01_annots_baseline_table.tsv", sep="\t", index_col=0)
tb2 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_annots_baseline_table.tsv", sep="\t", index_col=0)
tb3 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_annots_baseline_table.tsv", sep="\t", index_col=0)
tb1.columns = np.arange(0.1, 1.1, 0.1)
tb2.columns = np.arange(0.1, 1.1, 0.1)
tb3.columns = np.arange(0.1, 1.1, 0.1)
tb1 = tb1.iloc[[5,6,4,3,2,1,0],:]
tb2 = tb2.iloc[[5,6,4,3,2,1,0],:]
tb3 = tb3.iloc[[5,6,4,3,2,1,0],:]
f1 = tb1/tb1.sum()
f2 = tb2/tb2.sum()
f3 = tb3/tb3.sum()
print (f1)
print (f2)
print (f3)
#Plot in a single plot:
c7 = ["yellow", "tab:green", "yellowgreen", "tab:blue", "lightskyblue", "navy", "tab:gray"]
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb1.columns-0.03, tb1.sum(axis=0), log=True, width=0.03, color="tab:pink", label="JCTF (our study)")
ax[0].bar(tb3.columns, tb3.sum(axis=0), log=True, width=0.03, color="tab:olive", label="Zhang et al., AFR")
ax[0].bar(tb2.columns+0.03, tb2.sum(axis=0), log=True, width=0.03, color="tab:cyan", label="Zhang et al., EUR")
ax[0].set_ylabel("N(variant-gene)", fontsize=11)
ax[0].set_xlim([0.04,1.06])
ax[0].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Study:", fontsize=14)
ax[1].bar(f1.columns-0.03, f1.iloc[0,:], label=f1.index[0], color=c7[0], width=0.03)
for i in range(1,f1.shape[0]):
    ax[1].bar(f1.columns-0.03, f1.iloc[i,:], bottom=f1.iloc[:i,:].sum(axis=0), label=f1.index[i], color=c7[i], width=0.03)
ax[1].bar(f3.columns, f3.iloc[0, :], color=c7[0], width=0.03)
for i in range(1, f3.shape[0]):
    ax[1].bar(f3.columns, f3.iloc[i, :], bottom=f3.iloc[:i, :].sum(axis=0), color=c7[i], width=0.03)
ax[1].bar(f2.columns+0.03, f2.iloc[0,:], color=c7[0], width=0.03)
for i in range(1,f2.shape[0]):
    ax[1].bar(f2.columns+0.03, f2.iloc[i,:], bottom=f2.iloc[:i,:].sum(axis=0), color=c7[i], width=0.03)
#for ax[1], also the edge:
ax[1].bar(f1.columns-0.03, [1]*f1.shape[1], fill = False, edgecolor = "tab:pink", linewidth=2, width=0.028)
ax[1].bar(f3.columns+0.03, [1]*f3.shape[1], fill = False, edgecolor = "tab:cyan", linewidth=2, width=0.028)
ax[1].bar(f2.columns, [1]*f2.shape[1], fill = False, edgecolor = "tab:olive", linewidth=2, width=0.028)
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Variant class:", fontsize=14)
ax[1].set_xlabel("x$\leq$PIP", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticks(np.round(tb1.columns,2))
ax[1].set_xticklabels(np.round(tb1.columns,2))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_vep_vs_aa_and_ea_baseline.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aa_and_ea_baseline.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#EMS:
tb1 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pj_pip01_ems_table.tsv", sep="\t", index_col=0)
tb2 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/ea_pip01_ems_table.tsv", sep="\t", index_col=0)
tb3 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/aa_pip01_ems_table.tsv", sep="\t", index_col=0)
tb1.columns = np.arange(0.1, 1.1, 0.1)
tb2.columns = np.arange(0.1, 1.1, 0.1)
tb3.columns = np.arange(0.1, 1.1, 0.1)
f1 = tb1/tb1.sum()
f2 = tb2/tb2.sum()
f3 = tb3/tb3.sum()
#Plot in a single plot:
c6 = []
import matplotlib.cm as cm
vir = cm.viridis
for i in range(6):
    c6.append(vir(int(i/(5)*256)))
f1.index = ["(0.01,0.1]","(0.1,1]","(1,10]","(10,100]","(100,1000]","1000<"]
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(12, 6), gridspec_kw={'height_ratios': [1, 3]})
ax[0].bar(tb1.columns-0.03, tb1.sum(axis=0), log=True, width=0.03, color="tab:pink", label="JCTF (our study)")
ax[0].bar(tb3.columns, tb3.sum(axis=0), log=True, width=0.03, color="tab:olive", label="Zhang et al., AFR")
ax[0].bar(tb2.columns+0.03, tb2.sum(axis=0), log=True, width=0.03, color="tab:cyan", label="Zhang et al., EUR")
ax[0].set_ylabel("N(variant-gene)", fontsize=11)
ax[0].set_xlim([0.04,1.06])
ax[0].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Study:", fontsize=14)
ax[1].bar(f1.columns-0.03, f1.iloc[0,:], label=f1.index[0], color=c6[0], width=0.03)
for i in range(1,f1.shape[0]):
    ax[1].bar(f1.columns-0.03, f1.iloc[i,:], bottom=f1.iloc[:i,:].sum(axis=0), label=f1.index[i], color=c6[i], width=0.03)
ax[1].bar(f3.columns, f3.iloc[0, :], color=c6[0], width=0.03)
for i in range(1, f3.shape[0]):
    ax[1].bar(f3.columns, f3.iloc[i, :], bottom=f3.iloc[:i, :].sum(axis=0), color=c6[i], width=0.03)
ax[1].bar(f2.columns+0.03, f2.iloc[0,:], color=c6[0], width=0.03)
for i in range(1,f2.shape[0]):
    ax[1].bar(f2.columns+0.03, f2.iloc[i,:], bottom=f2.iloc[:i,:].sum(axis=0), color=c6[i], width=0.03)
#for ax[1], also the edge:
ax[1].bar(f1.columns-0.03, [1]*f1.shape[1], fill = False, edgecolor = "tab:pink", linewidth=2, width=0.028)
ax[1].bar(f3.columns+0.03, [1]*f3.shape[1], fill = False, edgecolor = "tab:cyan", linewidth=2, width=0.028)
ax[1].bar(f2.columns, [1]*f2.shape[1], fill = False, edgecolor = "tab:olive", linewidth=2, width=0.028)
ax[1].legend(bbox_to_anchor=(1.01,0.5), loc='center left', title="Expression modifier score (EMS)", fontsize=14)
ax[1].set_xlabel("x$\leq$PIP", fontsize=15)
ax[1].set_ylabel("Fraction", fontsize=15)
ax[1].set_ylim([0,1])
plt.subplots_adjust(hspace = .05)
ax[1].set_xticks(np.round(tb1.columns,2))
ax[1].set_xticklabels(np.round(tb1.columns,2))
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_vep_vs_aa_and_ea_ems.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_pip_vs_aa_and_ea_ems.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Vep annotations:
tb0 = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/vepannot_pqtl_baseline_stats_updated.tsv".format(chr), sep='\t', index_col=0)
annot_to_use = ["NA", "intron_variant", "upstream_gene_variant", "downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant",
              "synonymous_variant", "missense_variant","inframe", "lof"]
ys = {}
yerrs = {}
for annot in annot_to_use:
    tb = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/pqtl_enr_{0}_updated.tsv".format(annot), sep='\t', index_col=0)
    f0 = tb0.loc[True,annot] / tb0.loc[:,annot].sum()
    #fix the 0 bin
    toadd = list(np.array(tb0.loc[:,annot]) - np.array(tb.iloc[1:,:2].sum()))
    toadd.append(toadd[1]/(toadd[0]+toadd[1]))#adding the frac
    toadd.append(np.sqrt(toadd[2]*(1-toadd[2])/(toadd[0]+toadd[1])))
    tb.loc[0] = toadd
    tb.sort_index(inplace=True)
    ys[annot] = tb.frac/f0
    yerrs[annot] = tb.err / f0
ys = pd.DataFrame(ys)
yerrs = pd.DataFrame(yerrs)

c7 = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tk = ["[0,0.001)","[0.001,0.01)", "[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1]", "1"]
fig, ax = plt.subplots(1,3,sharex=False, sharey=True, figsize=(14,3))
idx = []
#negative guys:
ys.columns = ys.columns.str.replace("NA", "No annotations")
yerrs.columns = yerrs.columns.str.replace("NA", "No annotations")
ys_use = ys.loc[:,["No annotations", "intron_variant"]].T
yerrs_use = yerrs.loc[:,["No annotations", "intron_variant"]].T
x = np.arange(ys_use.shape[0])
i = 0
for pip in ys_use.columns:
    y = ys_use[pip]
    yerr = yerrs_use[pip]
    ax[0].errorbar(x-0.3+0.1*i, y, yerr, fmt="o", label=tk[i], color=c7[i])
    i += 1
idx = ["no annotation", "intron"]
ax[0].set_xticks(x, idx, rotation=30)
#nc guys
ys_use = ys.loc[:,["upstream_gene_variant", "downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"]].T
yerrs_use = yerrs.loc[:,["upstream_gene_variant", "downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant"]].T
x = np.arange(ys_use.shape[0])
i = 0
for pip in ys_use.columns:
    y = ys_use[pip]
    yerr = yerrs_use[pip]
    ax[1].errorbar(x-0.3+0.1*i, y, yerr, fmt="o", label=tk[i], color=c7[i])
    i += 1
idx = ["upstream gene\nvariant", "downstream gene\n variant", "3'UTR variant", "5'UTR variant"]
ax[1].set_xticks(x, idx, rotation=30)
#coding guys:
ys_use = ys.loc[:,["synonymous_variant", "missense_variant","inframe", "lof"]].T
yerrs_use = yerrs.loc[:,["synonymous_variant", "missense_variant","inframe", "lof"]].T
x = np.arange(ys_use.shape[0])
i = 0
for pip in ys_use.columns:
    y = ys_use[pip]
    yerr = yerrs_use[pip]
    ax[2].errorbar(x-0.3+0.1*i, y, yerr, fmt="o", label=tk[i], color=c7[i])
    i += 1
idx = ["synonymous\nvariant", "missense\n variant", "inframe", "LoF"]
ax[2].set_xticks(x, idx, rotation=30)
ax[0].set_xlabel("Annotation\n(Negative)")
ax[1].set_xlabel("Annotation\n(Non-coding)")
ax[2].set_xlabel("Annotation\n(Coding)")
ax[0].set_ylabel("Enrichment")
ax[0].legend(title="pQTL PIP bin")
ax[0].spines.right.set_visible(False)
ax[0].spines.top.set_visible(False)
ax[1].spines.right.set_visible(False)
ax[1].spines.top.set_visible(False)
ax[2].spines.right.set_visible(False)
ax[2].spines.top.set_visible(False)
ax[0].set_xlim([-0.4,3.4])
ax[1].set_xlim([-0.4,3.4])
ax[2].set_xlim([-0.4,3.4])
plt.subplots_adjust(wspace=0.05)
plt.yscale('log')
#ax[0].set_ylim([10**-1.1,10**4.1])
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_vep_enr_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_vep_enr_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()


#Uniprot annotations:
tbback = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/tmp/vepannot_background_uniprot_enr_updated.tsv", sep='\t', index_col=0)
tbback.index = [True, False, "frac", "err"] #str => bool
pj = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/vepannot_pip0001_aminoacids_uniprotannot_updated.tsv.gz", sep='\t')
interest = ["chains", "structure_beta", "structure_helix", "disulfbonds", "extracellular", "transmembrane", 'cytoplasmic']
tbs = {}
for categ in interest:
    f0 = tbback.loc["frac", categ]
    tb = pj.groupby(["min_pip_bin", categ]).size().unstack().fillna(0).astype(int)
    if tb.shape[1]==1:
        tb[True] = 0
    tb.loc[0] = tbback.iloc[[1, 0], :][categ] - tb.sum()
    tb.sort_index(inplace=True)
    tb["frac"] = tb[True] / tb.sum(axis=1)
    tb["err"] = np.sqrt(tb.frac*(1-tb.frac)/(tb[True]+tb[False]))
    tb["enr"] = tb.frac / f0
    tb["enr_err"] = tb.err / f0
    tbs[categ] = tb
    print (tb)
ys = {}
yerrs = {}
for annot in ["chains", "structure_beta", "structure_helix",  "extracellular", 'cytoplasmic']:
    tb = tbs[annot]
    ys[annot] = tb.enr
    yerrs[annot] = tb.enr_err
ys = pd.DataFrame(ys).sort_index()
yerrs = pd.DataFrame(yerrs).sort_index()

c7 = ["tab:gray", "tab:blue", "tab:green", "tab:olive", "tab:orange", "tab:red", "magenta"]
tk = ["[0,0.001)","[0.001,0.01)", "[0.01,0.1)", "[0.1,0.5)","[0.5,0.9)","[0.9,1)", "1"]
fig, ax = plt.subplots(1,2,sharex=False, sharey=True, figsize=(11,3))
idx = []
#structure guys:
ys_use = ys.loc[:,["chains", "structure_beta", "structure_helix"]].T
yerrs_use = yerrs.loc[:,["chains", "structure_beta", "structure_helix"]].T
x = np.arange(ys_use.shape[0])
i = 0
for pip in ys_use.columns:
    y = ys_use[pip]
    yerr = yerrs_use[pip]
    ax[0].errorbar(x-0.15+0.1*i, y, yerr, fmt="o", label=tk[i], color=c7[i])
    i += 1
idx = ["Chains", "beta-sheet", "Helix"]
ax[0].set_xticks(x, idx, rotation=30)
#location guys:
ys_use = ys.loc[:,["extracellular", 'cytoplasmic']].T
yerrs_use = yerrs.loc[:,["extracellular", 'cytoplasmic']].T
x = np.arange(ys_use.shape[0])
i = 0
for pip in ys_use.columns:
    y = ys_use[pip]
    yerr = yerrs_use[pip]
    ax[1].errorbar(x-0.15+0.1*i, y, yerr, fmt="o", label=tk[i], color=c7[i])
    #add arrow for too low/high
    for j in range(len(y)):
        if y[j] - yerr[j] < 1/10:
            ax[1].annotate('', xy=(x[j]-0.15+0.1*i, 1/10-0.00015), xytext=(x[j]-0.15+0.1*i, y[j]),
                           arrowprops=dict(arrowstyle='->', color=c7[i]))
    i += 1
idx = ["Extracellular", "Cytoplasmic"]
ax[1].set_xticks(x, idx, rotation=30)
ax[0].set_xlabel("Annotation\n(Structure)")
ax[1].set_xlabel("Annotation\n(Location)")
ax[0].set_ylabel("Enrichment")
ax[1].legend(title="pQTL PIP bin\n(Missense)")
ax[0].spines.right.set_visible(False)
ax[0].spines.top.set_visible(False)
ax[1].spines.right.set_visible(False)
ax[1].spines.top.set_visible(False)
plt.ylim([0.1,10])
ax[0].set_xlim([-0.4,2.7])
ax[1].set_xlim([-0.4,2.7])
plt.subplots_adjust(wspace=0.01)
plt.yscale("log")
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_strc_enr_updated2.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/pqtl_strc_enr_updated2.pdf', bbox_inches='tight', dpi=500)
plt.clf()

##eQTL-pQTL colocalization:

#Numbers:
st = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_table.tsv", sep='\t', index_col=0)
st.columns = ["[0.01,0.1)", "[0.1,0.9)", "[0.9,1]"]
st["total"] = st.sum(axis=1)
gtexcol = pd.read_json("~/Desktop/resources/gtex_colors.json").T
tissues = gtexcol.index.str.replace(" - ","_").str.replace(" \\(","_").str.replace("\\)","").str.replace(" ","_")
tissues = tissues.str.replace("Cells_Transformed_fibroblasts", "Cells_Cultured_fibroblasts") #manual transformation for this one alone
gtexcol.index = tissues
st = st.join(gtexcol, how="left")
st.sort_values(by="total", ascending=False, inplace=True)
colors = gtexcol.loc[st.index,"tissue_color_hex"]
colors = "#"+colors
#replace the brain colors with darker yellow: tab:olive is fine
colors = colors.str.replace("#EEEE00","tab:olive")
del st["total"]
st.index = st.index.str.replace("_"," ")
c3 = ["tab:brown", "tab:orange", "tab:red"]
ax = st.plot.bar(stacked=True, color=c3, figsize=(9,6))#, hatch=h3)
ax.legend(title='CLPP bin')
ax.set_xlabel("GTEx tissue")
ax.set_ylabel("Number of variant-gene pairs")
for i in range(st.shape[0]):
    ax.get_xticklabels()[i].set_color(colors[i])
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/coloc_gtextissues_num_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/coloc_gtextissues_num_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Corr. with causal var sample size
ncausal = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_n_pip09_pertissue.tsv", sep='\t', index_col=0)
sumclpp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_sum_pertissue.tsv", sep='\t', index_col=0)
ncausal.sort_index(inplace=True)
sumclpp.sort_index(inplace=True)
colors = gtexcol.loc[st.index,"tissue_color_hex"]
colors = "#"+colors
colors.sort_index(inplace=True)
print (sum(ncausal.index!=sumclpp.index), sum(ncausal.index!=colors.index))#making sure that the indices are the same
plt.figure(figsize=(4,3.3))
plt.scatter(ncausal, sumclpp, c=colors)
plt.xlabel("Number of variant-gene with GTEx PIP>0.9")
plt.ylabel("Sum of CLPP with plasma pQTL")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/coloc_gtextissues_num_scatter_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/coloc_gtextissues_num_scatter_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#Corr with the enrichment score is not that high:
enr = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_enrichment_scores.tsv", sep='\t', index_col=0)
enr.sort_index(inplace=True)
print (sum(ncausal.index!=enr.index))
plt.figure(figsize=(4,3.3))
plt.scatter(ncausal, enr.enr, c=colors)
plt.xlabel("Number of variant-gene with GTEx PIP>0.9")
plt.ylabel("Colocalization enrichment")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/coloc_gtextissues_enr_scatter_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/coloc_gtextissues_enr_scatter_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#And finally the coloc enrichment score plot:
enr = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/basic_stats/gtex_clpp_enrichment_scores.tsv", sep='\t', index_col=0)
gtexcol = pd.read_json("~/Desktop/resources/gtex_colors.json").T
tissues = gtexcol.index.str.replace(" - ","_").str.replace(" \\(","_").str.replace("\\)","").str.replace(" ","_")
tissues = tissues.str.replace("Cells_Transformed_fibroblasts", "Cells_Cultured_fibroblasts") #manual transformation for this one alone
gtexcol.index = tissues
enr = enr.join(gtexcol, how="left")
enr.sort_values(by="enr", inplace=True, ascending=False)
enr["idx"] = enr.tissue_abbrv
enr.loc[enr.pval<0.05,"idx"] = "*" + enr.loc[enr.pval<0.05,"idx"]
enr.loc[enr.pval<0.05/enr.shape[0],"idx"] = "*" + enr.loc[enr.pval<0.05,"idx"]
plt.figure(figsize=(8,3.5))
for x in range(enr.shape[0]): #Loop for different colors, although not ideal..
    y = enr.enr[x]
    yerr_upper = enr.upper[x] - y
    yerr_lower = y - enr.lower[x]
    col = "#" + enr.tissue_color_hex[x]
    plt.errorbar(x, y, [[yerr_lower], [yerr_upper]], color=col, fmt="o")
plt.axhline(y=1, linestyle=":", linewidth=0.5, zorder=-1.0, color="black")
plt.xticks(range(enr.shape[0]), enr.idx, rotation=90)
plt.xlabel("GTEx tissue")
plt.ylabel("Colocalization enrichment")
plt.xlim([-0.5,48.5])
plt.ylim([1/10,10])
plt.yscale('log')
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/gtex_coloc_enr_updated.png", bbox_inches='tight', dpi=500)
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/gtex_coloc_enr_updated.pdf", bbox_inches='tight', dpi=500)
plt.clf()




