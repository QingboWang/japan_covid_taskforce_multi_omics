#Plots used for constructing NPX matrix:
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


#1. PCA showing that no outliers are in the batch.
from sklearn.decomposition import PCA
df = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/n1404_parsed_matrix_invnorm.tsv.gz", sep='\t', compression='gzip', index_col=0, header=[0,1])
X = df.T
pca = PCA(n_components=10)
pca.fit(X)
print(pca.explained_variance_ratio_)
pcs = pd.DataFrame(pca.fit_transform(X))
pcs.index = X.index
vr = pca.explained_variance_ratio_
fig, ax = plt.subplots(1, 3, figsize=(12,3.8))
for i in range(3):
    ax[i].scatter(pcs[pcs.index.get_level_values(1)=="Q_00984"].iloc[:,i], pcs[pcs.index.get_level_values(1)=="Q_00984"].iloc[:,i+1], color="tab:green", label="Batch 1")
    ax[i].scatter(pcs[pcs.index.get_level_values(1)=="Q_01085"].iloc[:,i], pcs[pcs.index.get_level_values(1)=="Q_01085"].iloc[:,i+1], color="tab:red", label="Batch 2")
    ax[i].scatter(pcs[pcs.index.get_level_values(1)=="Q_01086"].iloc[:,i], pcs[pcs.index.get_level_values(1)=="Q_01086"].iloc[:,i+1], color="tab:blue", label="Batch 3")
    #ax[i].set_aspect('equal', adjustable='datalim')
    ax[i].set_xlabel("PC{0} ({1}%)".format(i+1, np.round(vr[i]*100, 3)))
    ax[i].set_ylabel("PC{0} ({1}%)".format(i+2, np.round(vr[i+1]*100, 3)))
ax[2].legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
fig.suptitle("Principal components (PCs) for protein expression matrix",y=0.92)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/npx_pca.png', bbox_inches='tight', dpi=500)
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/npx_pca.pdf', bbox_inches='tight', dpi=500)
plt.clf()

#2. Showing that concatenation of duplication works.
mat = pd.read_csv("~/Desktop/taskforce_n1102/olink_data/n1300_normed_npx.tsv", sep='\t')
#filter to CT cases
ids = pd.read_excel("~/Desktop/pqtl/Q-01085_Keio_2021/shirokane upload COVID-19 Taskforce manifest(Q-00984,Q-1085)_OsakaHV_modified_2022.04.20.xlsx")
ids_Q01086 = pd.read_excel("~/Desktop/taskforce_n1102/olink_data/Q-01086_Sample_Manifest.xlsx")
mat.index = mat.PlateID.str.replace("-","_").str.split("_").str[1]+"_"+mat.SampleID
ids.index = ids["Quote Number"].str.replace("-","_")+"_"+ids["Unique Sample ID \n(for letters standard A-Z only)"].fillna("NA").astype(str)
ids_Q01086.index = ids_Q01086["Quote Number"].str.replace("-","_")+"_"+ids_Q01086["Unique Sample ID \n(for letters standard A-Z only)"].fillna("NA").astype(str)
#雑に全てunderbarにこてい:
mat.index = "Q_"+mat.index.str.replace("-","_")
ids.index = ids.index.str.replace("-","_")
ids_Q01086.index = ids_Q01086.index.str.replace("-","_")
ctids = pd.concat([ids['Variable 2 (Subject)'],ids_Q01086['Variable 2 (Subject)']]).str.replace("-","_")
ctids.name = "CT_ID"
mat = mat.join(ctids, how="left")
mat = mat[mat.CT_ID.fillna("NA").str.startswith("CT_")]
dup_proteins = mat.Assay.value_counts()[mat.Assay.value_counts()>1420]
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats
sns.set(font_scale=1.5)
for gene in ['CXCL8', 'IDO1', 'IL6', 'LMOD1', 'TNF', 'SCRIB']:
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

#3. Showing that duplicated samples have okay replication rate:
mat85 = mat[mat.index=="Q_01085_NUPK131384"].sort_values(by="OlinkID")
mat86 = mat[mat.index=="Q_01086_NUPK131384"].sort_values(by="OlinkID")
mat85.index = mat85.OlinkID
mat86.index = mat86.OlinkID
mat85 = mat85.join(mat86, how="left", rsuffix="86")
vc = mat.CT_ID.value_counts()
dupid = vc[vc==vc.max()].index
i = 0
j = 0
cnt = 1
fix, ax = fig, ax = plt.subplots(4, 4, figsize=(12, 12))
for id in dupid:
    mat85 = mat[(mat.CT_ID==id)&(~mat.index.str.startswith("Q_01086"))].sort_values(by="OlinkID")
    mat86 = mat[(mat.CT_ID==id)&(mat.index.str.startswith("Q_01086"))].sort_values(by="OlinkID")
    mat85.index = mat85.OlinkID
    mat86.index = mat86.OlinkID
    mat85 = mat85.join(mat86, how="left", rsuffix="86")
    x = mat85.NPX
    y = mat85.NPX86
    ax[i,j].scatter(x, y, alpha=0.4)
    ax[i,j].axvline(x=0, linestyle="--", linewidth=0.5, color="black")
    ax[i,j].axhline(y=0, linestyle="--", linewidth=0.5, color="black")
    m = mat85.NPX.min()
    M = mat85.NPX.max()
    ax[i,j].plot([m, M], [m, M], linestyle="--", linewidth=0.5, color="black")
    ax[i,j].set_title("Sample ID: {0}, r={1}".format(cnt, np.round(stats.pearsonr(x[(~x.isna())&(~y.isna())], y[(~x.isna())&(~y.isna())])[0], 3)), fontsize=12)
    if i==3:
        ax[i,j].set_xlabel("Batch 1 or 2")
    if j==0:
        ax[i,j].set_ylabel("Batch 3")
    if i<3:
        i = i + 1
    else:
        j = j + 1
        i = 0
    cnt += 1
fig.suptitle("Correlation between duplicated samples in different batch", y=0.95)
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig("/Users/qingbowang/Desktop/plots/dupsamples_normalization_check.png", dpi=400)
plt.savefig("/Users/qingbowang/Desktop/plots/dupsamples_normalization_check.pdf", dpi=400)
plt.clf()


