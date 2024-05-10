import pandas as pd
import numpy as np

#NPX files
ol1 = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/Q-00984_NPX.csv", sep=';')
ol2 = pd.read_csv("~/Desktop/pqtl/Q-01085_Keio_2021/Q-01085_Namkoong_EDTAPlasma_NPX_2022-03-18.csv", sep=';') #start from un-normalized guy?
ol3 = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/Q-01086_Namkoong_NPX_2022-09-22.csv", sep=';')

#Names (CT ID)
ol12_names = pd.read_excel("~/Desktop/pqtl/Q-01085_Keio_2021/shirokane upload COVID-19 Taskforce manifest(Q-00984,Q-1085)_OsakaHV_modified_2022.04.20.xlsx")
ol12_names.index = ol12_names["Unique Sample ID \n(for letters standard A-Z only)"]
ol3_names = pd.read_excel("~/Desktop/taskforce_n1102/olink_data/COVID-19 Taskforce JPN-SM_20220818.xlsx")
ol3_names.index = ol3_names["Unique Sample ID \n(for letters standard A-Z only)"]

#join the CT_ID
print ("before join: {0}, {1}, {2}".format(ol1.shape, ol2.shape, ol3.shape))
ol1.index = ol1.SampleID
ol1 = ol1.join(ol12_names[["Variable 2 (Subject)"]], how="left") #OL1 somehow already has those
ol2.index = ol2.SampleID
ol2 = ol2.join(ol12_names[["Variable 2 (Subject)","Remarks"]], how="left")
ol3.index = ol3.SampleID
ol3 = ol3.join(ol3_names[["Variable 2 (Subject)","Remarks"]], how="left")
print ("after join: {0}, {1}, {2}".format(ol1.shape, ol2.shape, ol3.shape))
ol1["Variable 2 (Subject)"].unique()
ol2["Variable 2 (Subject)"].unique()
ol3["Variable 2 (Subject)"].unique()
ol1[ol1["Variable 2 (Subject)"].isna()].SampleID.unique()
ol2[ol2["Variable 2 (Subject)"].isna()].SampleID.unique()
ol3[ol3["Variable 2 (Subject)"].isna()].SampleID.unique()

#Check that for OL1, no discrepancy between the original sample ID and that we annotated:
ol1[ol1.SubjectID != ol1["Variable 2 (Subject)"]][["SubjectID", "Variable 2 (Subject)"]].drop_duplicates(keep="first") #NA anyways, so this is good.

#change to lower
ol1["Variable 2 (Subject)"] = ol1["Variable 2 (Subject)"].str.replace("-","_")
ol2["Variable 2 (Subject)"] = ol2["Variable 2 (Subject)"].str.replace("-","_")
ol3["Variable 2 (Subject)"] = ol3["Variable 2 (Subject)"].str.replace("-","_")

#check the remarks:
ol1.Remarks.unique()
ol2.Remarks.unique()
ol3.Remarks.unique()

#and save this updated csv, with same sep=";"
n1 = len(ol1.SampleID.unique())
n2 = len(ol2.SampleID.unique())
n3 = len(ol3.SampleID.unique())
print ("before filtering non-CT samples: {0}, {1}, {2}".format(n1, n2, n3))
ol1.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-00984_NPX_idannot.csv", sep=';', index=False)
ol2.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-01085_NPX_idannot.csv", sep=';', index=False)
ol3.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-01086_NPX_idannot.csv", sep=';', index=False)

#and also the ones where non-CT samples are removed
ol1filt = ol1[ol1["Variable 2 (Subject)"].fillna("NA").str.startswith("CT_")]
ol2filt = ol2[ol2["Variable 2 (Subject)"].fillna("NA").str.startswith("CT_")]
ol3filt = ol3[ol3["Variable 2 (Subject)"].fillna("NA").str.startswith("CT_")]
n1 = len(ol1filt.SampleID.unique())
n2 = len(ol2filt.SampleID.unique())
n3 = len(ol3filt.SampleID.unique())
print ("After filtering non-CT samples: {0}, {1}, {2}".format(n1, n2, n3)) #1ズレているのはNAが消えただけと思われる..
ol1filt.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-00984_NPX_idannot_filtered_to_CT.csv", sep=';', index=False)
ol2filt.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-01085_NPX_idannot_filtered_to_CT.csv", sep=';', index=False)
ol3filt.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-01086_NPX_idannot_filtered_to_CT.csv", sep=';', index=False)

ol1filt.groupby(["Variable 2 (Subject)", "Assay_Warning"]).size() #Assay warning is per protein.
ol1filt.groupby(["Variable 2 (Subject)", "QC_Warning"]).size()
ol1filt[ol1filt.Assay_Warning!="PASS"].NPX.sort_values()
ol1filt[ol1filt.QC_Warning!="PASS"].NPX.sort_values()

#Fill these with NA actually
ol1filt2 = ol1filt.copy(deep=True)
ol1filt2.loc[(ol1filt2.Assay_Warning!="PASS")|(ol1filt2.QC_Warning!="PASS"),"NPX"] = np.nan
ol2filt2 = ol2filt.copy(deep=True)
ol2filt2.loc[(ol2filt2.Assay_Warning!="PASS")|(ol2filt2.QC_Warning!="PASS"),"NPX"] = np.nan
ol3filt2 = ol3filt.copy(deep=True)
ol3filt2.loc[(ol3filt2.Assay_Warning!="PASS")|(ol3filt2.QC_Warning!="PASS"),"NPX"] = np.nan

ol1filt2.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-00984_NPX_warningfilt_filtered_to_CT.csv", sep=';', index=False)
ol2filt2.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-01085_NPX_warningfilt_filtered_to_CT.csv", sep=';', index=False)
ol3filt2.to_csv("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-01086_NPX_warningfilt_filtered_to_CT.csv", sep=';', index=False)

#And plug these simply to R for olink normalization; done
mat = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/n1300_normed_npx_20230114.tsv", sep='\t')

#samples:
mat["Variable.2..Subject."].value_counts() #1413 now

#Additionally mask blacklist genes for Q01086
blacklist_genes = ["ALMS1","ARAF","BRD1","CD82",
    "CEACAM20","CEP290","COL4A4","CPLX2",
    "CTAG1A_CTAG1B","EVPL","FUT1","GABARAPL1",
    "IFIT1","ITGAX","MAGEA3","MTHFSD",
    "OGT","RAPGEF2","SAP18","SEPTIN7",
    "SH3GL3","TAGLN3","TAP1","TET2",
    "TNPO1","TPPP2","UPK3BL1","ZNF75D"]
mat.loc[mat.Assay.apply(lambda x: x in blacklist_genes) & mat.PlateID.str.startswith("Q_01086"),'NPX'] = np.nan #takes some time

#First, filter by fraction of NAs per gene:
gnafrac = mat.groupby("Assay").NPX.apply(lambda x: sum(x.isna()/len(x))).sort_values()
gnafrac.tail(20) #easy, 25% should be a good threshold
genes_to_filt = gnafrac[gnafrac>0.25].index

#fraction of NAs per sample:
nafrac = mat.groupby("Variable.2..Subject.").NPX.apply(lambda x: sum(x.isna()/len(x))).sort_values()
nafrac.tail(50)
nafrac.hist(bins=20)
plt.show()

#get the list of redundant genes:
mat.Assay.value_counts().value_counts()
dup_proteins = mat.Assay.value_counts()[mat.Assay.value_counts()==5716].index
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats
for gene in dup_proteins:
    sub = mat[(mat.Assay==gene)]
    sub.sort_values(by=["CT_ID", "Panel"], inplace=True)
    sub.reset_index(inplace=True)
    df = pd.concat([sub.iloc[::4,:].NPX.reset_index(drop=True), sub.iloc[1::4,:].NPX.reset_index(drop=True),
                    sub.iloc[2::4,:].NPX.reset_index(drop=True), sub.iloc[3::4,:].NPX.reset_index(drop=True)], axis=1)
    df.columns = sub.Panel.head(4) #hard-coding but fine...
    rs = pd.Series(df.columns).apply(lambda x: df.apply(lambda y: stats.pearsonr(df[x][(~df[x].isna()&(~y.isna()))], y[(~df[x].isna()&(~y.isna()))])[0]))
    minr = rs.min().min()
    sns.pairplot(df)
    plt.suptitle("NPX in 4 assays for gene {0}\nr$\geq${1}".format(gene, minr))
    plt.tight_layout()
    plt.savefig("/Users/qingbowang/Desktop/plots/npx_duplication_{0}.png".format(gene), dpi=500)
    plt.savefig("/Users/qingbowang/Desktop/plots/npx_duplication_{0}.pdf".format(gene), dpi=500)
    plt.clf()
vc = mat.CT_ID.value_counts()
dupid = vc[vc==vc.max()].index
i = 0
j = 0
fix, ax = fig, ax = plt.subplots(4, 4, figsize=(12, 12))
for id in dupid:
    mat85 = mat[(~mat.PlateID.str.startswith("Q_01086")) & (mat.CT_ID == id)].sort_values(by="OlinkID")
    mat86 = mat[mat.PlateID.str.startswith("Q_01086") & (mat.CT_ID == id)].sort_values(by="OlinkID")
    mat85.index = mat85.OlinkID
    mat86.index = mat86.OlinkID
    mat85 = mat85.join(mat86, how="left", rsuffix="86")
    mat85 = mat85[~(mat85.NPX.isna()|mat85.NPX86.isna())]
    x = mat85.NPX
    y = mat85.NPX86
    ax[i,j].scatter(x, y, alpha=0.4)
    ax[i,j].axvline(x=0, linestyle="--", linewidth=0.5, color="black")
    ax[i,j].axhline(y=0, linestyle="--", linewidth=0.5, color="black")
    m = mat85.NPX.min()
    M = mat85.NPX.max()
    ax[i,j].plot([m, M], [m, M], linestyle="--", linewidth=0.5, color="black")
    ax[i,j].set_title("{0}, r={1}".format(id, np.round(stats.pearsonr(x[(~x.isna())&(~y.isna())], y[(~x.isna())&(~y.isna())])[0], 3)), fontsize=11)
    if i==3:
        ax[i,j].set_xlabel("Q_00984 or Q_01085")
    if j==0:
        ax[i,j].set_ylabel("Q_01086")
    if i<3:
        i = i + 1
    else:
        j = j + 1
        i = 0
plt.savefig("/Users/qingbowang/Desktop/plots/dupsamples_normalization_check.png", dpi=500)
plt.savefig("/Users/qingbowang/Desktop/plots/dupsamples_normalization_check.pdf", dpi=500)
plt.clf()
#filter the 86 one is fine.
dupid = vc[vc==vc.max()].index
matfilt = mat[~(mat.PlateID.str.startswith("Q_01086")&mat.CT_ID.apply(lambda x: x in dupid))]

#make it a sample-gene matrix at this stage:
#to do so, we would like to keep the plateID:
matfilt["batch"] = "NA" #dummy
matfilt.loc[matfilt.PlateID.str.startswith("Q-00984"),"batch"] = "Q_00984"
matfilt.loc[matfilt.PlateID.str.startswith("Q-01085"),"batch"] = "Q_01085"
matfilt.loc[matfilt.PlateID.str.startswith("Q_01086"),"batch"] = "Q_01086"
matfilt.batch.value_counts()
#and now we can make the matrix:
tb = matfilt.groupby(["CT_ID", "batch", "Assay"]).NPX.median().unstack()
tb.to_csv("~/Desktop/taskforce_n1102/olink_data/n1413_parsed_matrix.tsv", sep='\t') #the full matrix is saved (invnorm not done yet)
#Remember; at this stage, no filtering based on nafrac is done yet.

#and from here we first need to filter based on relatedness status:
rel = pd.read_csv("~/Desktop/taskforce_n1102/n1300/relatedness_from_chr1_n1419.tsv.relatedness2", sep="\t", index_col=[0,1])
rl = rel.RELATEDNESS_PHI.unstack()
keep = np.triu(np.ones(rl.shape)).astype('bool').reshape(rl.size)
rl = rl.stack()[keep]
rl = rl[rl.index.get_level_values(0)!=rl.index.get_level_values(1)] #1004653 = 1418*1417/2 yes.
idt = rl[rl>0.48] #identical pair
sib = rl[(rl>0.15)&(rl<0.48)] #making lower - let this be the sibling pairs

#first, filter out the ones where genotype do not exist
genotyped_samples = rel.index.get_level_values(0).unique()
tb["batch"] = tb.index.get_level_values(1)
tb = tb.droplevel(1)
tb = tb.loc[tb.index.intersection(genotyped_samples),:]
tb.set_index("batch", append=True, inplace=True) #1401 genotyped samples pre-QC

#then remove the identical samples based on NA rates:
nafrac = tb.isna().sum(axis=1)/tb.shape[1]
nafrac.index = nafrac.index.get_level_values(0)
nafrac = pd.DataFrame(nafrac)
nafrac.columns = ["na_frac"]
nafrac.index.names = ["INDV1"]
#annotate nafrac for identical samples
idt = pd.DataFrame(idt)
idt = idt.join(nafrac.na_frac, how="left", on="INDV1")
nafrac.index.names = ["INDV2"]
idt = idt.join(nafrac.na_frac, how="left", on="INDV2", rsuffix="_indv2")
#annotate nafrac for siblings
sib = pd.DataFrame(sib)
nafrac.index.names = ["INDV1"]
sib = sib.join(nafrac.na_frac, how="left", on="INDV1")
nafrac.index.names = ["INDV2"]
sib = sib.join(nafrac.na_frac, how="left", on="INDV2", rsuffix="_indv2")
idt = idt.fillna(1) #1 = when the sample do not exist
#to_remove = samples to remove for NAfrac being relatively large
to_remove_idt = idt[idt.na_frac>idt.na_frac_indv2].index.get_level_values(0).union(idt[idt.na_frac<idt.na_frac_indv2].index.get_level_values(1))
to_remove_sib = sib[sib.na_frac>sib.na_frac_indv2].index.get_level_values(0).union(sib[sib.na_frac<sib.na_frac_indv2].index.get_level_values(1))
#remove these siblings:
tb["batch"] = tb.index.get_level_values(1)
tb = tb.droplevel(1)
tb = tb.loc[tb.index.difference(to_remove_idt.union(to_remove_sib)),:]
tb.set_index("batch", append=True, inplace=True) #1401 genotyped samples pre-QC

#and now we do the NAs filtering:
nafrac = (tb.isna().sum(axis=1)/tb.shape[1])
nafrac_gene = (tb.isna().sum(axis=0)/tb.shape[0])
nafrac.sort_values().tail(100)#Let's make 33.333% as the threshold (arbitral...)
nafrac_gene.sort_values().tail(100) #Let's make 33.333% as the threshold
tb = tb.loc[nafrac<1/3,nafrac_gene<1/3] #Here we go, 1369 samples passing stringent QC. #No we actually went on with 1384 samples, filling NAs
tb.to_csv("~/Desktop/taskforce_n1102/olink_data/n1369_postfilt_parsed_matrix.tsv.gz", sep='\t', compression='gzip') #note: still keeping NAs
#inverse normal transformation:
import pandas as pd
from scipy import stats
def inverse_normal_transform(M):
    """Transform rows to a standard normal distribution"""
    if isinstance(M, pd.Series):
        r = stats.rankdata(M)
        return pd.Series(stats.norm.ppf(r/(M.shape[0]+1)), index=M.index, name=M.name)
    else:
        R = stats.rankdata(M, axis=1)  # ties are averaged
        Q = stats.norm.ppf(R/(M.shape[1]+1))
        if isinstance(M, pd.DataFrame):
            Q = pd.DataFrame(Q, index=M.index, columns=M.columns)
        return Q
tb = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/n1369_postfilt_parsed_matrix.tsv.gz", sep='\t', index_col=[0,1]).T #gene x sample
tb_meanimp = tb.fillna(tb.mean(axis=1))
tb_meanimp_invnorm = inverse_normal_transform(tb_meanimp)
tb_meanimp_invnorm = tb_meanimp_invnorm.T
for genepairs in tb_meanimp_invnorm.loc[:,tb_meanimp_invnorm.columns.str.contains("_")]:
    for gene in genepairs.split("_"):
        tb_meanimp_invnorm.loc[:,gene] = tb_meanimp_invnorm.loc[:,genepairs]
    del tb_meanimp_invnorm[genepairs]
tb_meanimp_invnorm.T.to_csv("~/Desktop/taskforce_n1102/olink_data/n1369_postfilt_parsed_matrix_invnorm.tsv.gz", sep='\t', compression='gzip')

#PCA:
tb = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/n1369_postfilt_parsed_matrix_invnorm.tsv.gz", sep='\t', compression='gzip',
                 index_col=0, header=[0,1]).T
from sklearn.decomposition import PCA
X = tb
pca = PCA(n_components=10)
pca.fit(X)
print(pca.explained_variance_ratio_)
pcs = pd.DataFrame(pca.fit_transform(X))
pcs.index = X.index
pcs.to_csv("~/Desktop/taskforce_n1102/olink_data/n1369_pcs.tsv", sep='\t')

plt.scatter(pcs[pcs.index.get_level_values(1)=="Q_00984"].iloc[:,0], pcs[pcs.index.get_level_values(1)=="Q_00984"].iloc[:,1], color="tab:green")
plt.scatter(pcs[pcs.index.get_level_values(1)=="Q_01085"].iloc[:,0], pcs[pcs.index.get_level_values(1)=="Q_01085"].iloc[:,1], color="tab:red")
plt.scatter(pcs[pcs.index.get_level_values(1)=="Q_01086"].iloc[:,0], pcs[pcs.index.get_level_values(1)=="Q_01086"].iloc[:,1], color="tab:blue")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.show()

plt.scatter(pcs[pcs.index.get_level_values(1)=="Q_00984"].iloc[:,2], pcs[pcs.index.get_level_values(1)=="Q_00984"].iloc[:,3], color="tab:green")
plt.scatter(pcs[pcs.index.get_level_values(1)=="Q_01085"].iloc[:,2], pcs[pcs.index.get_level_values(1)=="Q_01085"].iloc[:,3], color="tab:red")
plt.scatter(pcs[pcs.index.get_level_values(1)=="Q_01086"].iloc[:,2], pcs[pcs.index.get_level_values(1)=="Q_01086"].iloc[:,3], color="tab:blue")
plt.xlabel("PC3")
plt.ylabel("PC4")
plt.show()

#also check the D_i (after filtering NAs)
df = tb.T
corrmat = pd.DataFrame(index=df.columns, columns = df.columns)
for i in range(corrmat.shape[0]):
    l = df.iloc[:, i]
    corrmat.iloc[:,i] = df.apply(lambda x: stats.pearsonr(x, l)[0], axis=0)
    if i%10==0:
        print ("done {0}, {1}".format(i, tm.ctime())) #takes a few 5~10 minutes.
corrmat.to_csv("~/Desktop/taskforce_n1102/olink_data/n1369_corrmat.tsv", sep='\t')
r_i = (corrmat.sum(axis=0)-1)/(corrmat.shape[1]-1)
r_barbar = (corrmat.sum().sum()-corrmat.shape[0])/(corrmat.shape[0]**2-corrmat.shape[0])
D_i = (r_i - r_barbar)/(np.median(r_i - r_barbar))
D_i.to_csv("~/Desktop/taskforce_n1102/olink_data/n1369_corr_D_i.tsv", sep='\t')
#plot:
D_i = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/n1369_corr_D_i.tsv", sep='\t', index_col=[0,1], header=[0,1], squeeze=True)
mu = D_i.mean()
sigma = D_i.std()
plt.figure(figsize=(4.5,3))
p = D_i.hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    if ((rectangle.get_x() < mu - 2 * sigma)): #High -> OK
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=mu-2*sigma, linestyle="--", linewidth=1, color="black")
plt.axvline(x=mu+2*sigma, linestyle="--", linewidth=1, color="black")
plt.axvline(x=mu, linestyle=":", linewidth=1, color="black")
plt.xlabel("Intersample correlation deviation ($D_i$)")
plt.ylabel("Number of samples")
plt.show()

#add start and end etc from gencode (considering the strand of course):
mat = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/n1369_postfilt_parsed_matrix_invnorm.tsv.gz", sep='\t', compression='gzip',index_col=0, header=[0,1])
ensgids = pd.read_csv("~/Desktop/taskforce_n1102/gencode.v30.genes.parsed.tsv", sep='\t')
ensgids.index = ensgids.gene_name
mat.index.difference(ensgids.index) #Those that are in different name in gencode v30
#manually annotated non-matching ones:
nonmatch = pd.read_csv("~/Desktop/taskforce_n1102/manualfilled_nonmatch.tsv", sep='\t', header=None)
nonmatch.columns = ["pname", "rname"]
nonmatch.rname = nonmatch.rname.str.replace("_","") #calendar対策
nonmatch.rname = nonmatch.rname.str.replace(" ","") #その他誤字対策
nonmatch.loc[nonmatch.rname=="NPPB","rname"] = "??"
#remove the ones with ??
nonmatch = nonmatch[nonmatch.rname!="??"]
nonmatch.index = nonmatch.pname
#change the gene_name in pqtl mat accordingly:
def conv_func(x):
    try: return (nonmatch.loc[x, "rname"])
    except: return (x)
mat.index = pd.Series(mat.index).apply(lambda x: conv_func(x))
still_missing = np.setdiff1d(mat.index, ensgids.index) #those can be removed
mat = mat.loc[np.setdiff1d(mat.index, still_missing),:]

meta = ensgids.loc[mat.index, ["#CHROM","start","end","ensg_id", "strand"]]
#shxt there are some non-unique element... also some manual cleaning here:
meta = meta[meta.ensg_id!="ENSG00000284934.1"]
meta['idx'] = meta.index
meta.drop_duplicates(subset="idx", inplace=True)
del meta["idx"]
#now what are the guys missing in meta?? No it is just about duplicates
np.setdiff1d(mat.index, meta.index) #Now everything is uniquely matched.
mat.columns = mat.columns.get_level_values(0)
out = meta.join(mat, how="inner")
out.rename(columns = {"#CHROM":"#chr","ensg_id":"gene_id"}, inplace=True)
out["start"] = out.start-1
out.loc[out.strand=="-","start"] = out.loc[out.strand=="-","end"]-1
out["end"] = out.start + 1
del out["strand"]
#Now the start and end is either (TSS+1, TES+1) or (TES+1, TSS+1) depending on strand
out.sort_values(by=["#chr", "start"], inplace=True)
out.to_csv("~/Desktop/taskforce_n1102/n1369.protein.expression.bed", sep='\t', index=False)
out.to_csv("~/Desktop/taskforce_n1102/n1369.protein.expression.indexon.bed", sep='\t', index=True)

#JCTF eQTL expression matrix as a reference (and to make sure that the start and end are identical)
eq = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.bed.gz", sep="\t")
eq = eq.iloc[:,:4]
pq = pd.read_csv("~/Desktop/taskforce_n1102/n1369.protein.expression.bed", sep='\t')
pq = pq.iloc[:,:4]
eq.index = eq.gene_id
pq.index = pq.gene_id
intsct = eq.index.intersection(pq.index)
eq.loc[intsct,:]
pq.loc[intsct,:]
sum(eq.loc[intsct,"start"]!=pq.loc[intsct,"start"]) #is zero. Perfect.

# prep. covariates

#genotype PCs
#in local
df = pd.read_csv("~/Downloads/n1384.protein.expression.bed", sep="\t")
df = pd.DataFrame(df.columns[4:])
df[1] = 0
df = df.iloc[:,[1,0]]
df.to_csv("~/Downloads/n1384_protein_samplenames_w_famid.txt", sep="\t", header=None, index=False)
#in bash:
"""protein_pca.sh
cd n1419_pqtl_call
n1384list=/home/qwang/n1419_pqtl_call/n1384_protein_samplenames_w_famid.txt #create this
vcf=/home/qwang/n1419_vcf_hg38/taskforce_n1419_imputed.dose.mac3.rsq6.hg38.samplefilt.vcf.gz
../plink_dir/plink --pca --vcf $vcf --allow-extra-chr --const-fid 0 --keep $n1384list
"""
import pandas as pd
df = pd.read_csv("/home/qwang/n1419_pqtl_call/plink.eigenvec", sep=" ", index_col=1, header=None).T
df = df.iloc[1:,:] #removing the first redundant line
df = df.iloc[:5,:] #5 PCs is enough
df.index = ["PC1","PC2","PC3","PC4","PC5"]
df.columns.name = "ID"
df.to_csv("/home/qwang/n1419_pqtl_call/n1384_forprotein_genotype_pcs.tsv", sep="\t")

#Other covariates such as severity are omitted in this code due to privacy