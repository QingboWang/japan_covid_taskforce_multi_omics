import pandas as pd
import numpy as np
import time as tm
from matplotlib import pyplot as plt
from scipy import stats
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib.colors import LogNorm
import seaborn as sns


#get the lead variant based on marginal stats:
p0 = []
e0 = []
for chr in list(range(1,23))+["X"]:
    p = pd.read_csv("~/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t', compression="gzip")
    e = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t', compression="gzip")
    p0.append(p.sort_values(by=["gene_id", "pval_nominal"], ascending=True).groupby("gene_id").head(1))
    e0.append(e.sort_values(by=["gene_id", "pval_nominal"], ascending=True).groupby("gene_id").head(1))
    print ("done{0}, {1}".format(chr, tm.ctime()))
p0 = pd.concat(p0)
e0 = pd.concat(e0)
p0.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/trans_qtl_revision/cis_pqtl_lead_vars.tsv", sep="\t")
e0.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/trans_qtl_revision/cis_eqtl_lead_vars.tsv", sep="\t")

pos = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/trans_qtl_revision/cis_eqtl_lead_vars.tsv", sep="\t")
pos.index = pos.variant_id
pos.index.name = "ID"
pos["chrom"] = pos.variant_id.str.split(":").str[0]
pos["pos"] = pos.variant_id.str.split(":").str[1].astype(int)
pos = pos.iloc[:,-2:]
pos.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/for_trans_qtl_lead_eqtls.variant_df.vcf.gz", sep="\t", compression="gzip")

pos = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/trans_qtl_revision/cis_pqtl_lead_vars.tsv", sep="\t")
pos.index = pos.variant_id
pos.index.name = "ID"
pos["chrom"] = pos.variant_id.str.split(":").str[0]
pos["pos"] = pos.variant_id.str.split(":").str[1].astype(int)
pos = pos.iloc[:,-2:]
pos.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/for_trans_qtl_lead_pqtls.variant_df.vcf.gz", sep="\t", compression="gzip")

#And get the genotypes as well:
#(already done for neg)
import pandas as pd
import numpy as np
import time as tm
#To select the raws:
pose = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/for_trans_qtl_lead_eqtls.variant_df.vcf.gz", sep="\t", compression="gzip", index_col=0)
posp = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/for_trans_qtl_lead_pqtls.variant_df.vcf.gz", sep="\t", compression="gzip", index_col=0)
#To select the cols:
st = pd.read_csv("~/Desktop/taskforce_n1102/n1369.protein.expression.bed", sep='\t') #already inv-normed and is n1384
sample_protein = st.columns[4:]
st = pd.read_csv("~/Desktop/taskforce_n1102/n1019.expression.bed.gz", sep='\t') #Actually not protein but fine to name prot here, for simplicity
sample_rna = st.columns[4:]
import io
import gzip
def read_vcf_zipped(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
dfpos_e = []
dfpos_p = []
for chr in list(range(1,23))+["X"]:
    df = read_vcf_zipped("/Users/qingbowang/Desktop/taskforce_n1102/n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_chr{0}_reIDed.vcf.gz".format(chr))
    df.index = df.ID
    dfpos_e_sub = df.loc[df.index.intersection(pose.index),:]
    dfpos_p_sub = df.loc[df.index.intersection(posp.index), :]
    dfpos_e.append(dfpos_e_sub.loc[:,sample_rna])
    dfpos_p.append(dfpos_p_sub.loc[:,sample_protein])
    print ("Done {0}, {1}".format(chr, tm.ctime())) #Takes time but will be done in a few hours
dfpos_e = pd.concat(dfpos_e)
dfpos_p = pd.concat(dfpos_p)
dfpos_e = dfpos_e.applymap(lambda x: float(x.split(":")[1]))
dfpos_p = dfpos_p.applymap(lambda x: float(x.split(":")[1]))
dfpos_e.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/trans_qtl_revision/for_trans_qtl_lead_eqtls.genotype_df.vcf.gz", sep="\t", compression="gzip")
dfpos_p.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/trans_qtl_revision/for_trans_qtl_lead_pqtls.genotype_df.vcf.gz", sep="\t", compression="gzip")




#and go for trans-QTL call:
#Done: mkdir /home/qwang/n1419_trans
#module use /usr/local/package/modulefiles/
#module load python/3.8
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans
#pQTL:
phenotype_bed_file = '/home/qwang/n1419_pqtl_call/n1384.protein.expression.bed.gz'
covariates_file = '/home/qwang/n1419_pqtl_call/n1384.protein.combined_covariates_idadded.txt'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T  # samples x covariates
#case:
fn = "/home/qwang/n1419_trans/for_trans_qtl_lead_pqtls.genotype_df.vcf.gz"
genotype_df = pd.read_csv(fn, sep="\t", index_col=0)
fn = "/home/qwang/n1419_trans/for_trans_qtl_lead_pqtls.variant_df.vcf.gz"
variant_df = pd.read_csv(fn, sep="\t", index_col=0)
prefix = "pqtl_lead" #gw for genome wide
trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, return_sparse=True, pval_threshold=1, maf_threshold=0.01, batch_size=20000)#i.e. no any pval threshold
trans_df.to_csv("/home/qwang/n1419_trans/n1384.trans.pqtl.pvals.lead.tsv.gz", sep="\t", compression="gzip")
trans_df = trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=5000000)
trans_df.to_csv("/home/qwang/n1419_trans/n1384.trans.pqtl.pvals.lead.cisfiltered.tsv.gz", sep="\t", compression="gzip")
#cont - already done:

#eQTL:
phenotype_bed_file = '/home/qwang/n602_plus_500_parseddata/n1019.expression.bed.gz'
covariates_file = "/home/qwang/n602_plus_500_parseddata/n1019.sevon.combined_covariates_idadded.txt"
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T  # samples x covariates
#case:
fn = "/home/qwang/n1419_trans/for_trans_qtl_lead_eqtls.genotype_df.vcf.gz"
genotype_df = pd.read_csv(fn, sep="\t", index_col=0)
fn = "/home/qwang/n1419_trans/for_trans_qtl_lead_eqtls.variant_df.vcf.gz"
variant_df = pd.read_csv(fn, sep="\t", index_col=0)
prefix = "eqtl_lead" #gw for genome wide
trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df,
                           return_sparse=True, pval_threshold=1e-2, maf_threshold=0.01,
                           batch_size=20000)
trans_df.to_csv("/home/qwang/n1419_trans/n1019.trans.eqtl.pvals.lead.tsv.gz", sep="\t", compression="gzip")
trans_df = trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=5000000)
trans_df.to_csv("/home/qwang/n1419_trans/n1019.trans.eqtl.pvals.lead.cisfiltered.tsv.gz", sep="\t", compression="gzip")
#cont - already done

#For plot: move things to local.
#QQ Plot suggesting that cis-pQTLs are good candidate of trans-pQTLs:
case = pd.read_csv("~/Downloads/n1384.trans.pqtl.pvals.lead.cisfiltered.tsv.gz", sep="\t", index_col=0)
cont = pd.read_csv("~/Downloads/n1384.trans.pqtl.pvals.cont.cisfiltered.tsv.gz", sep="\t", index_col=0)
case.sort_values(by="pval", inplace=True, ascending=True)
case["expected_pval"] = np.arange(1,case.shape[0]+1)/case.shape[0]
case["x"] = -np.log10(case.expected_pval)
case["y"] = -np.log10(case.pval)

cont.sort_values(by="pval", inplace=True, ascending=True)
cont["expected_pval"] = np.arange(1,cont.shape[0]+1)/cont.shape[0]
cont["x"] = -np.log10(cont.expected_pval)
cont["y"] = -np.log10(cont.pval)
from scipy.stats import chi2
l_case = chi2.ppf(1-np.median(case.pval), 1)/chi2.ppf(0.5, 1)
l_cont = chi2.ppf(1-np.median(cont.pval), 1)/chi2.ppf(0.5, 1)
plt.figure(figsize=(4,5))
plt.scatter(cont.x, cont.y, color="tab:blue", label="Random\n($\lambda={0}$)".format(np.round(l_cont,4)))
plt.scatter(case.x, case.y, color="tab:orange", label="Lead cis-pQTL\n($\lambda={0}$)".format(np.round(l_case,4)))
plt.plot([0,10],[0,10], linestyle="--", linewidth=1, color="black")
plt.xlim([0,10])
plt.xlabel("Expected trans-pQTL -log10(p)")
plt.ylabel("Observed trans-pQTL -log10(p)")
plt.legend(title="Variants:", loc="upper left")
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/trans_pqtl_qq.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/trans_pqtl_qq.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#overview plot
df = pd.read_csv("~/Downloads/n1384.trans.pqtl.pvals.lead.tsv.gz", sep="\t", index_col=0)
df.set_index("phenotype_id", inplace=True)
#add the position of the trans-gene (TSS):
genes = pd.read_csv('~/Downloads/n1384.protein.expression.bed.gz', sep="\t", usecols=[0,1,2,3])
genes.set_index("gene_id", inplace=True)
genes["gene_center"] = (genes.start + genes.end) / 2
genes.columns = ["gene_chr", "start", "end", "gene_center"]
df = df.join(genes[["gene_chr", "start", "gene_center"]], how="left")

#adding cis-gene name:
cisg = pd.read_csv("~/Downloads/cis_pqtl_lead_vars.tsv", sep="\t")
cisg.set_index("variant_id", inplace=True)
df.index.name = "gene_id"
df.reset_index(inplace=True)
df.set_index("variant_id", inplace=True)
df = df.join(cisg.gene_id, how="left", rsuffix="_cis")
df.reset_index(inplace=True)
#creating the absolute position:
import pysam
ref = pysam.FastaFile('/Users/qingbowang/Downloads/hg38.fa')
chromosomes = ref.references
chromosome_lengths = {}
for chromosome in chromosomes:
    length = ref.get_reference_length(chromosome)
    chromosome_lengths[chromosome] = length
chrs = [str(i) for i in range(1, 23)] + ["X"]
len_to_add = {}
len_to_add[chrs[0]] = 0
for i in range(1,len(chrs)):
    len_to_add[chrs[i]] = len_to_add[chrs[i-1]] + chromosome_lengths["chr{0}".format(chrs[i-1])]
xlab = []
for i in range(len(chrs)):
    xlab.append(len_to_add[chrs[i]] + chromosome_lengths["chr{0}".format(chrs[i])]/2) #=center of each chr
df["variant_chr"] = df.variant_id.str.split(":").str[0]
df["variant_pos"] = df.variant_id.str.split(":").str[1].astype(int)
df["variant_pos_toadd"] = df.variant_chr.str.replace("chr","").apply(lambda x: len_to_add[x])
df["variant_abs_pos"] = df.variant_pos + df.variant_pos_toadd
df["gene_pos_toadd"] = df.gene_chr.str.replace("chr","").apply(lambda x: len_to_add[x])
df["gene_abs_pos"] = df.gene_center + df.gene_pos_toadd
df["gene_abs_tss"] = df.start + df.gene_pos_toadd

n_test_trans = df[abs(df.variant_abs_pos-df.gene_abs_tss)>5*10**6].shape[0]

#adding gene name:
gns = pd.read_csv("~/Downloads/gene_names.txt", sep="\t", index_col=0)
df.index = df.gene_id.str.split("\\.").str[0]
df = df.join(gns["Gene name"], how="left")
df.index = df.gene_id_cis.str.split("\\.").str[0]
df = df.join(gns["Gene name"], how="left", rsuffix=" cis")
all_trans = df[abs(df.variant_abs_pos-df.gene_abs_tss)>5*10**6]
all_signif = all_trans[all_trans.pval<0.05/all_trans.shape[0]]
cisgenes = all_signif["Gene name cis"].value_counts() #Happy with ABO result.
cisgenes_lenient = all_trans[all_trans.pval<1e-5]["Gene name cis"].value_counts()
to_focus = cisgenes[cisgenes>1]
all_signif.index = all_signif["Gene name cis"].fillna("NA")
to_focus_pos = all_signif.loc[to_focus.index,:].drop_duplicates(subset="Gene name cis").variant_abs_pos
to_focus = pd.concat([to_focus, to_focus_pos], axis=1)
#Manual cocatenation:
cisg = pd.read_csv("~/Downloads/cis_pqtl_lead_vars.tsv", sep="\t")
cisg.index = cisg.gene_id.str.split("\\.").str[0]
cisg = cisg.join(gns, how="left")
cisg.index = cisg["Gene name"]
cisg.loc[["ABO","OBP2B","ITIH1", "ITIH4"],"pval_nominal"]

to_focus.columns = ["N", "pos"]

import matplotlib.cm as cm
vir = cm.viridis
dfplt = df[df.pval < 1e-5]
dfcis = dfplt[abs(dfplt.variant_abs_pos-dfplt.gene_abs_tss)<5*10**6]
dftrans = dfplt[abs(dfplt.variant_abs_pos-dfplt.gene_abs_tss)>5*10**6]
plt.figure(figsize=(12,8))
plt.scatter(dfcis.variant_abs_pos, dfcis.gene_abs_pos, color="tab:grey", edgecolor="black", s=20, zorder=3, clip_on=False, label="Cis")
tr = dftrans[dftrans.pval>0.05/n_test_trans] #1e-8 -ish
plt.scatter(tr.variant_abs_pos, tr.gene_abs_pos, color=vir(0), edgecolor="black", s=20, zorder=3, clip_on=False, label="p<1e-5")
tr = dftrans[(dftrans.pval<0.05/n_test_trans)&(dftrans.pval>1e-20)]
plt.scatter(tr.variant_abs_pos, tr.gene_abs_pos, color=vir(0.4), edgecolor="black", s=35, zorder=3, clip_on=False, label="p<1.1e-8\n(Bonferroni)")
tr = dftrans[(dftrans.pval<1e-20)&(dftrans.pval>1e-40)]
plt.scatter(tr.variant_abs_pos, tr.gene_abs_pos, color=vir(0.65), edgecolor="black", s=50, zorder=3, clip_on=False, label="p<1e-20")
tr = dftrans[(dftrans.pval<1e-40)]
plt.scatter(tr.variant_abs_pos, tr.gene_abs_pos, color=vir(0.9), edgecolor="black", s=65, zorder=3, clip_on=False, label="p<1e-40")
for pos in len_to_add.values():
    plt.axvline(x=pos, linestyle="--", linewidth=0.5, color="black", zorder=-2)
    plt.axhline(y=pos, linestyle="--", linewidth=0.5, color="black", zorder=-2)
#adding gene annots:
for i in range(to_focus.shape[0]):
    x = to_focus.pos.values[i]
    txtx = x
    s = to_focus.index.values[i] + " (n={0})".format(to_focus.N.values[i])
    if to_focus.index.values[i]=="ABO":
        txtx = x - (10**8)/3
    elif to_focus.index.values[i]=="OBP2B":
        txtx = x + (10**8)/4
    elif to_focus.index.values[i]=="ITIH1":
        txtx = x+(10**8)/2
    plt.axvline(x=x, linestyle="-", linewidth=0.5, color="tab:green", zorder=-2)
    plt.text(x=txtx, y=len_to_add["X"] + chromosome_lengths["chrX"], s=s, rotation=60, va="bottom", ha="left")
plt.xticks(xlab, list(range(1,23))+["X"], rotation=45)
plt.yticks(xlab, list(range(1,16))+["","","","...","","",""]+["X"], rotation=45)
plt.xlabel("Variant position (Chromosome)", fontsize=14)
plt.ylabel("Gene position (Chromosome)", fontsize=14)
plt.xlim([len_to_add["1"], len_to_add["X"]+chromosome_lengths["chrX"]])
plt.ylim([len_to_add["1"], len_to_add["X"]+chromosome_lengths["chrX"]])
plt.legend(loc='upper left', bbox_to_anchor=(0, 1.4))
plt.tight_layout()
plt.savefig('/Users/qingbowang/Downloads/trans_pqtl_overview.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Downloads/trans_pqtl_overview.pdf', bbox_inches='tight', dpi=500)
plt.clf()


#HLA focused:

gns = pd.read_csv("~/Desktop/resources/gene_names.txt", sep="\t", index_col=0)
#trans eQTL
trans_e = pd.read_csv("~/Desktop/taskforce_n1102/n1005.hla.eqtl.pvals.1e2.cisfiltered.tsv", sep="\t") #broad capturing
trans_e.index = trans_e.phenotype_id.str.split("\\.").str[0]
trans_e = trans_e.join(gns["Gene name"], how="left")
trans_e.sort_values(by="pval").head(50) #KIRとかは共通してそうだが..
trans_e.set_index(["variant_id", "phenotype_id"], inplace=True)
trans_e_fourd = trans_e[trans_e.index.get_level_values(0).str.count(":")==1]
trans_e_fourd["hla_gene"] = np.array(pd.Series(trans_e_fourd.index.get_level_values(0).str.split("_").str[:-1]).apply(lambda x: "_".join(x)))

#trans pQTL
trans_df = pd.read_csv("~/Desktop/taskforce_n1102/n1005.hla.pvals.1e2.cisfiltered.tsv", sep="\t") #broad capturing
trans_df.index = trans_df.phenotype_id.str.split("\\.").str[0]
trans_df = trans_df.join(gns["Gene name"], how="left")
trans_df.set_index(["variant_id", "phenotype_id"], inplace=True)
trans_p_fourd = trans_df[trans_df.index.get_level_values(0).str.count(":")==1]
trans_p_fourd["hla_gene"] = np.array(pd.Series(trans_p_fourd.index.get_level_values(0).str.split("_").str[:-1]).apply(lambda x: "_".join(x)))

#list of tested genes:
mat = pd.read_csv("~/Desktop/taskforce_n1102/n1300/mrna_and_protein_n998.tsv.gz", sep="\t", compression="gzip", index_col=0)
genes = mat.columns[:int(len(mat.columns)/2)] #half are just "_protein"
genes_unv = pd.Series(genes).str.split("\\.").str[0]
measured_gene_names = gns.loc[gns.index.intersection(genes_unv),"Gene name"]

#e.g. all:
cols = {}
cols["HLA-A"] = "tab:orange"
cols["HLA-B"] = "tab:red"
cols["HLA-C"] = "tab:pink"
cols["HLA-DRB1"] = "slateblue"
cols["HLA-DQA1"] = "tab:cyan"
cols["HLA-DQB1"] = "tab:blue"
cols["HLA-DPA1"] = "lightgreen"
cols["HLA-DPB1"] = "tab:green"
cols["MICA"] = "tab:olive"
shapes = {}
shapes["HLA-A"] = "o"
shapes["HLA-B"] = "o"
shapes["HLA-C"] = "o"
shapes["HLA-DRB1"] = "D"
shapes["HLA-DQA1"] = "D"
shapes["HLA-DQB1"] = "D"
shapes["HLA-DPA1"] = "D"
shapes["HLA-DPB1"] = "D"
shapes["MICA"] = "s"
hla_genes = ["HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "MICA"] #These are the genes of interest after all

#get the position:
genes_bed = pd.read_csv('/Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.bed.gz', sep="\t")
genes_bed = genes_bed.iloc[:,:4]
genes_bed.index = genes_bed.gene_id
genes_bed.index.names = ["phenotype_id"]
trans_e_fourd = trans_e_fourd.join(genes_bed, how="left", on="phenotype_id")
genes_bed = pd.read_csv('/Users/qingbowang/Desktop/taskforce_n1102/n1005.protein.expression.tsscorrect.bed', sep="\t")
genes_bed = genes_bed.iloc[:,:4]
genes_bed.index = genes_bed.gene_id
genes_bed.index.names = ["phenotype_id"]
trans_p_fourd = trans_p_fourd.join(genes_bed, how="left", on="phenotype_id")
right_edge = {} #right edge for each chr
for chr in list(range(1,23))+["X"]:
    chr = str(chr)
    right_edge[chr] = max(trans_e_fourd[trans_e_fourd["#chr"]=="chr"+chr].start.max(), trans_p_fourd[trans_p_fourd["#chr"]=="chr"+chr].start.max())
cum_size = {}
chrs = list(range(1,23))+["X"]
cum_size["chr1"] = 0
for i in range(1,len(chrs)):
    cum_size["chr"+str(chrs[i])] = cum_size["chr"+str(chrs[i-1])] + right_edge[str(chrs[i-1])]
cum_size["NA"] = 0


trans_e_fourd["cum_pos"] = trans_e_fourd["#chr"].apply(lambda x: cum_size[x])
trans_e_fourd["abs_pos"] = trans_e_fourd.start + trans_e_fourd.cum_pos
trans_p_fourd["cum_pos"] = trans_p_fourd["#chr"].fillna("NA").apply(lambda x: cum_size[x])
trans_p_fourd["abs_pos"] = trans_p_fourd.start + trans_p_fourd.cum_pos

trans_e_fourd["hla_gene"] = trans_e_fourd["hla_gene"].str.replace("_","-")
trans_p_fourd["hla_gene"] = trans_p_fourd["hla_gene"].str.replace("_","-")

trans_p_fourd = trans_p_fourd.loc[trans_p_fourd["#chr"]!="chrX",:]
trans_e_fourd = trans_e_fourd.loc[trans_e_fourd["#chr"]!="chrX",:]

toplot_p = trans_p_fourd.sort_values(by="pval").groupby(["hla_gene","Gene name"]).head(1) #各 hla class - trans geneに関して最小値だけ残せば良い.
toplot_e = trans_e_fourd.sort_values(by="pval").groupby(["hla_gene","Gene name"]).head(1) #各 hla class - trans geneに関して最小値だけ残せば良い.
plt.rcParams.update({'font.size': 14})
xlab_pos = []
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(15, 5))
for i in range(22):
    if i%2==0: c='#bbbbbbff'
    else: c = '#eeeeeeff'
    ax[0].axvspan(cum_size["chr"+str(chrs[i])], cum_size["chr"+str(chrs[i+1])], facecolor=c, zorder=-2) #color for manhattan
    ax[1].axvspan(cum_size["chr" + str(chrs[i])], cum_size["chr" + str(chrs[i + 1])], facecolor=c, zorder=-2)  # color for manhattan
    xlab_pos.append((cum_size["chr"+str(chrs[i])]+cum_size["chr"+str(chrs[i+1])])/2)
#White out the cis region
left = 25726063
right = 33400644
ax[0].axvspan(cum_size["chr6"]+left, cum_size["chr6"]+right, facecolor="#ffffffff", hatch='////', edgecolor='black', linewidth=0.5, zorder=-2) #color for manhattan
ax[1].axvspan(cum_size["chr6"]+left, cum_size["chr6"]+right, facecolor="#ffffffff", hatch='////', edgecolor='black', linewidth=0.5, zorder=-2) #color for manhattan

ax[0].axhline(y=5, linestyle="--", linewidth=0.5, color="black", zorder=-1)
ax[1].axhline(y=5, linestyle="--", linewidth=0.5, color="black", zorder=-1)
#for g in hla_classes:
for g in hla_genes: #restricting to 1-ish guys only
    ax[0].scatter(toplot_p[toplot_p.hla_gene==g].abs_pos,
                  -np.log10(toplot_p[toplot_p.hla_gene==g].pval), color=cols[g], marker=shapes[g], s=6, edgecolors="#000000ff", linewidths=0.1)
    ax[1].scatter(toplot_e[toplot_e.hla_gene==g].abs_pos,
                  -np.log10(toplot_e[toplot_e.hla_gene==g].pval), color=cols[g], marker=shapes[g], s=6, edgecolors="#000000ff", linewidths=0.1)
    ax[0].scatter(toplot_p[(toplot_p.pval<1e-5)&(toplot_p.hla_gene==g)].abs_pos,
                  -np.log10(toplot_p[(toplot_p.pval<1e-5)&(toplot_p.hla_gene==g)].pval), color=cols[g], marker=shapes[g], label=g, s=24, edgecolors="#000000ff", linewidths=0.2)
    ax[1].scatter(toplot_e[(toplot_e.pval<1e-5)&(toplot_e.hla_gene==g)].abs_pos,
                  -np.log10(toplot_e[(toplot_e.pval<1e-5)&(toplot_e.hla_gene==g)].pval), color=cols[g], marker=shapes[g], s=24, edgecolors="#000000ff", linewidths=0.2)
#gene annotation:
gannot_p = toplot_p[(toplot_p.pval<1e-5)&(toplot_p.hla_gene.apply(lambda x: x in hla_genes))].groupby("Gene name").head(1)
gannot_e = toplot_e[(toplot_e.pval<1e-5)&(toplot_e.hla_gene.apply(lambda x: x in hla_genes))].groupby("Gene name").head(1)
gannot_p["coord_x"] = gannot_p.abs_pos
gannot_e["coord_x"] = gannot_e.abs_pos #based on this, change when there are too many duplicates
gannot_p.loc[gannot_p["Gene name"]=="GP1BB", "coord_x"] = gannot_p.loc[gannot_p["Gene name"]=="GP1BB", "coord_x"]-10**7.8
gannot_e.loc[gannot_e["Gene name"]=="SLCO4A1", "coord_x"] = gannot_e.loc[gannot_e["Gene name"]=="SLCO4A1", "coord_x"]-10**7.5
gannot_e.loc[gannot_e["Gene name"]=="CHRNA2", "coord_x"] = gannot_e.loc[gannot_e["Gene name"]=="CHRNA2", "coord_x"]+10**7
gannot_e.loc[gannot_e["Gene name"]=="NRG1", "coord_x"] = gannot_e.loc[gannot_e["Gene name"]=="NRG1", "coord_x"]+10**7
gannot_e.loc[gannot_e["Gene name"]=="TRBV27", "coord_x"] = gannot_e.loc[gannot_e["Gene name"]=="TRBV27", "coord_x"]-10**8.1
gannot_e.loc[gannot_e["Gene name"]=="TRBV9", "coord_x"] = gannot_e.loc[gannot_e["Gene name"]=="TRBV9", "coord_x"]-10**8.1
gannot_e.loc[gannot_e["Gene name"]=="TRBV15", "coord_x"] = gannot_e.loc[gannot_e["Gene name"]=="TRBV15", "coord_x"]-10**8.1
gannot_e.loc[gannot_e["Gene name"]=="TRBV3-1", "coord_x"] = gannot_e.loc[gannot_e["Gene name"]=="TRBV3-1", "coord_x"]-10**8.1
gannot_p["coord_y"] = -np.log10(gannot_p.pval)+1
gannot_e["coord_y"] = -np.log10(gannot_e.pval)+0.5 #based on this, change when there are too many duplicates
gannot_p.loc[gannot_p["Gene name"]=="ICAM2", "coord_y"] = gannot_p.loc[gannot_p["Gene name"]=="ICAM2", "coord_y"]-3
gannot_p.loc[gannot_p["Gene name"]=="NPTX1", "coord_y"] = gannot_p.loc[gannot_p["Gene name"]=="NPTX1", "coord_y"]+2.5
gannot_p.loc[gannot_p["Gene name"]=="GP1BB", "coord_y"] = gannot_p.loc[gannot_p["Gene name"]=="GP1BB", "coord_y"]+8
gannot_e.loc[gannot_e["Gene name"]=="RNF5P1", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="RNF5P1", "coord_y"]-8
gannot_e.loc[gannot_e["Gene name"]=="PFKFB3", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="PFKFB3", "coord_y"]+1.5
gannot_e.loc[gannot_e["Gene name"]=="ZNF239", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="ZNF239", "coord_y"]+0.75
gannot_e.loc[gannot_e["Gene name"]=="CLC", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="CLC", "coord_y"]+20
gannot_e.loc[gannot_e["Gene name"]=="ZNF235", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="ZNF235", "coord_y"]+15
gannot_e.loc[gannot_e["Gene name"]=="TTYH1", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="TTYH1", "coord_y"]+10
gannot_e.loc[gannot_e["Gene name"]=="KIR2DL1", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="KIR2DL1", "coord_y"]+4
gannot_e.loc[gannot_e["Gene name"]=="CHRNA2", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="CHRNA2", "coord_y"]+4
gannot_e.loc[gannot_e["Gene name"]=="TRBV27", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="TRBV27", "coord_y"]
gannot_e.loc[gannot_e["Gene name"]=="TRBV9", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="TRBV9", "coord_y"]+2.5
gannot_e.loc[gannot_e["Gene name"]=="TRBV15", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="TRBV15", "coord_y"]+5
gannot_e.loc[gannot_e["Gene name"]=="TRBV3-1", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="TRBV3-1", "coord_y"]+7.5
gannot_e.loc[gannot_e["Gene name"]=="PFKFB3", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="PFKFB3", "coord_y"]+1
gannot_e.loc[gannot_e["Gene name"]=="ZNF239", "coord_y"] = gannot_e.loc[gannot_e["Gene name"]=="ZNF239", "coord_y"]+1.1
for i in range(gannot_p.shape[0]):
    if gannot_p["Gene name"][i] in measured_gene_names.values:
        ax[0].text(gannot_p.coord_x[i], gannot_p.coord_y[i], gannot_p["Gene name"][i], rotation=45, fontsize=10, weight="bold")
    else:
        ax[0].text(gannot_p.coord_x[i], gannot_p.coord_y[i], gannot_p["Gene name"][i], rotation=45, fontsize=10)
for i in range(gannot_e.shape[0]):
    if gannot_e["Gene name"][i] in measured_gene_names.values:
        ax[1].text(gannot_e.coord_x[i], gannot_e.coord_y[i], gannot_e["Gene name"][i], rotation=-45, fontsize=10, ma="left", va="top", fontstyle="oblique", weight="bold")
    else:
        ax[1].text(gannot_e.coord_x[i], gannot_e.coord_y[i], gannot_e["Gene name"][i], rotation=-45, fontsize=10, ma="left", va="top", fontstyle="oblique")
ax[0].text(10**7, 30, "Protein", color="black", fontsize=16)
ax[1].text(10**7, 30, "mRNA", color="black", fontsize=16)
ax[0].text(cum_size["chr6"]+(left+right)/2, 38, "(MHC)", color="black", fontsize=10, horizontalalignment = 'center', verticalalignment = 'center') ##あとはこれの位置調整.. TO DO..
ax[0].set_ylim([0,36])
ax[1].set_ylim([0,36])
ax[1].invert_yaxis()
ax[0].set_ylabel("trans-pQTL -log$_{10}$(p)")
ax[1].set_ylabel("trans-eQTL -log$_{10}$(p)")
ax[1].set_xlabel("Chromosome")
ax[1].set_xticks(xlab_pos)
ax[1].set_xticklabels(list(np.arange(1,22+1)))
ax[0].set_xlim([-10**6.5,max(toplot_e.abs_pos.max(), toplot_p.abs_pos.max())+10**6.5])
ax[0].legend(bbox_to_anchor=(1.01,0.2), loc='center left', title="Gene:")
plt.tight_layout()
plt.subplots_adjust(hspace = .1)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/hla_miami_updated.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/hla_miami_updated.pdf', bbox_inches='tight', dpi=500)
plt.clf()
