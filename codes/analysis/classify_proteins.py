#classify proteins based on blood mRNA expression level (for n=2244 intersections)
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

#Classify proteins based on whole blood expression
df = pd.read_csv("~/Desktop/taskforce_n1102/n1300/mrna_and_protein_n998.tsv.gz", sep="\t", compression="gzip", index_col=0)
dfm = df.loc[:,~df.columns.str.contains("protein")]
dfp = df.loc[:,df.columns.str.contains("protein")]
dfp.columns = dfm.columns
gt = pd.read_csv("~/Desktop/resources/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
del gt["Description"]
gt.index = gt.Name.str.split("\\.").str[0]
del gt["Name"]
wb_frac = gt["Whole Blood"]/gt.sum(axis=1)
wb_frac = wb_frac[~wb_frac.index.duplicated(keep='first')]
wb_frac = wb_frac.reindex(dfm.columns.str.split("\\.").str[0],fill_value=0)
wb_frac.name = "wbfrac"
data = pd.DataFrame(wb_frac)
data["wbfrac_bin"] = -1
data.loc[data.wbfrac<0.01, "wbfrac_bin"] = 0
data.loc[data.wbfrac>=0.01,"wbfrac_bin"] = 1
data.loc[data.wbfrac>0.1,"wbfrac_bin"] = 10
data.loc[data.wbfrac>0.25,"wbfrac_bin"] = 25
data.loc[data.wbfrac>0.5,"wbfrac_bin"] = 50
data = data[data.wbfrac_bin!=-1]
data.wbfrac_bin.value_counts()

#Plot ridge:
corr_all = pd.Series(dfm.columns).apply(lambda x: stats.pearsonr(dfm.loc[:,x], dfp.loc[:,x])[0])
corr_all.index = dfm.columns
corr_all.index = corr_all.index.str.split("\\.").str[0]
corr_all = pd.DataFrame(corr_all)
corr_all.columns = ["pearsonr"]
data = data.join(corr_all, how="left")
low = data[data.wbfrac_bin==0]
#Save:
data.to_csv("~/Desktop/taskforce_n1102/n1300/wbfrac_classification.tsv.gz", sep="\t", compression="gzip")

#Box plots
wb = pd.read_csv("~/Desktop/taskforce_n1102/n1300/wbfrac_classification.tsv.gz", sep="\t", compression="gzip", index_col=0)
wb["bin_detail"] = -1
wb.loc[wb.wbfrac>=0,"bin_detail"] = 0
wb.loc[wb.wbfrac>=0.01,"bin_detail"] = 1
wb.loc[wb.wbfrac>=0.05,"bin_detail"] = 5
wb.loc[wb.wbfrac>=0.25,"bin_detail"] = 25
wb.loc[wb.wbfrac>=0.50,"bin_detail"] = 50
from scipy import stats
df = pd.read_csv("~/Desktop/taskforce_n1102/n1300/mrna_and_protein_n998.tsv.gz", sep="\t", compression="gzip", index_col=0)
dfm = df.loc[:,~df.columns.str.contains("protein")]
dfp = df.loc[:,df.columns.str.contains("protein")]
dfp.columns = dfp.columns.str.replace("_protein", "")
meta = pd.read_csv("~/Desktop/taskforce_n1102/n1300/eqtl_n998.combined_covariates_idadded.txt", sep="\t", index_col=0)
meta = meta.iloc[-3:,:].T
dfm_sev = pd.concat([dfm.loc[dfm.index.intersection(meta[meta.severity>3].index), :], dfm.loc[dfm.index.intersection(meta[meta.severity==3].index), :].sample(n=1, random_state=0)], axis=0)
dfp_sev = dfp.loc[dfm_sev.index, :]
dfm_mild = dfm.loc[dfm.index.intersection(meta[meta.severity<3].index), :]
dfp_mild = dfp.loc[dfp.index.intersection(meta[meta.severity<3].index), :]
corr_sev = pd.Series(dfm_sev.columns).apply(lambda x: stats.pearsonr(dfm_sev.loc[:,x], dfp_sev.loc[:,x])[0])
corr_mild = pd.Series(dfm_mild.columns).apply(lambda x: stats.pearsonr(dfm_mild.loc[:,x], dfp_mild.loc[:,x])[0])
corr_mild.index = dfm.columns
corr_sev.index = dfm.columns
corrs = pd.Series(dfm.columns).apply(lambda x: stats.pearsonr(dfm.loc[:,x], dfp.loc[:,x])[0])
corr_sev = pd.Series(dfm_sev.columns).apply(lambda x: stats.pearsonr(dfm_sev.loc[:,x], dfp_sev.loc[:,x])[0])
corr_mild = pd.Series(dfm_mild.columns).apply(lambda x: stats.pearsonr(dfm_mild.loc[:,x], dfp_mild.loc[:,x])[0])
corrs.index = dfm.columns.str.split("\\.").str[0]
corr_sev.index = dfm_sev.columns.str.split("\\.").str[0]
corr_mild.index = dfm_mild.columns.str.split("\\.").str[0]
corrs.name = "corr"
corr_sev.name = "corr_sev"
corr_mild.name = "corr_mild"
wb = wb.join(corrs, how="left").join(corr_sev, how="left").join(corr_mild, how="left")
wb.sort_values(by="corr")
data = wb.copy(deep=True)
data.sort_values(by='wbfrac', ascending=True, inplace=True)
# Calculate the rolling mean and SEM
rolling_mean_x = data['wbfrac'].rolling(window=50, min_periods=50, center=True).mean()
rolling_sem_x = data['wbfrac'].rolling(window=50, min_periods=50, center=True).sem()
rolling_mean_y = data['corr'].rolling(window=50, min_periods=50, center=True).mean()
rolling_sem_y = data['corr'].rolling(window=50, min_periods=50, center=True).sem()
palette = sns.cubehelix_palette(5, start=2)
fig, ax = plt.subplots(figsize=(5,3))
# Plot the rolling mean
ax.plot(np.log10(rolling_mean_x), rolling_mean_y, color="black")
ax.fill_between(np.log10(rolling_mean_x), rolling_mean_y - rolling_sem_y, rolling_mean_y + rolling_sem_y, color="black", edgecolor="none", alpha=0.4)
ax.set_xlabel("%Whole Blood expression in GTEx\n(rolling average)")
ax.set_ylabel("Corr(mRNA, protein) in\npearson R (rolling average)")
ax.set_xticks([-4, -3, -2, np.log10(0.05), np.log10(0.25), np.log10(0.5), 0])
ax.set_xticklabels(["0.01","0.1","1","5", "25", "50", "100"])
ax.axvline(x=np.log10(0.01), color="tab:gray", linestyle="--", linewidth=1)
ax.axvline(x=np.log10(0.05), color="tab:gray", linestyle="--", linewidth=1)
ax.axvline(x=np.log10(0.25), color="tab:gray", linestyle="--", linewidth=1)
ax.axvline(x=np.log10(0.50), color="tab:gray", linestyle="--", linewidth=1)
plt.axvspan(np.log10(0.00004),np.log10(0.01), alpha=0.5, color=palette[0], zorder=-2)
plt.axvspan(np.log10(0.01),np.log10(0.05), alpha=0.5, color=palette[1], zorder=-2)
plt.axvspan(np.log10(0.05),np.log10(0.25), alpha=0.5, color=palette[2], zorder=-2)
plt.axvspan(np.log10(0.25),np.log10(0.5), alpha=0.5, color=palette[3], zorder=-2)
plt.axvspan(np.log10(0.5),np.log10(1), alpha=0.5, color=palette[4], zorder=-2)
plt.tight_layout()
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/expression_corr_rolling_ave.png', dpi=500, bbox_inches='tight')
plt.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/expression_corr_rolling_ave.pdf', dpi=500, bbox_inches='tight')
plt.clf()

#Ridge plot:
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
# Define the data and variables
df = wb.reset_index().sort_values(by="bin_detail")
cnt = df.bin_detail.value_counts()
x_var = "corr"
y_var = "bin_detail"
means = df.groupby(y_var)[x_var].mean()
row_order = wb.bin_detail.sort_values().unique()
palette = sns.cubehelix_palette(5, start=2)
palette_dict = dict(zip(row_order, palette))
titles = ["<1", "[1,5)","[5,25)","[25,50)","50$\leq$"]
# Create the ridge plot
plt.figure(figsize=(12,8))
g = sns.FacetGrid(df, row=y_var, hue=y_var, aspect=10, height=.5, palette=palette_dict, row_order=row_order)
g.map(sns.kdeplot, x_var, clip_on=False, fill=True, alpha=1, lw=1.5, bw_method=.25, clip=[-0.41, 0.81])
g.map(sns.kdeplot, x_var, clip_on=False, color="w", lw=2, bw_method=.25, clip=[-0.41, 0.81])
g.map(plt.axhline, y=0, lw=2, clip_on=False)
g.map(lambda **kwargs: plt.axvline(0, linestyle='--', color='tab:gray'), alpha=0.5, lw=1)
for i, ax in enumerate(g.axes.flat):
    ax.axvline(x=means[row_order[i]], linestyle='--', color='black')
    ax.text(x=-0.35, y=2, s=titles[i], color=palette_dict[row_order[i]])
    ax.text(x=0.5, y=2, s="{0} genes".format(int(cnt[row_order[i]])), color=palette_dict[row_order[i]])
    sns.despine(top=True, ax=ax)  # Remove the subplot title
g.set(yticks=[], ylabel=None, xlabel="Corr(mRNA, protein) in Pearson R")
g.set_titles("")
g.axes[-1, 0].set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
g.axes[2, 0].set_ylabel("Density", zorder=10)
g.axes[0, 0].text(x=-0.40, y=4, s="%Whole blood\nexpression (GTEx)", fontsize=10)
g.axes[0, 0].text(x=0.45, y=4.8, s="Number of genes", fontsize=11)
plt.xlim([-0.41, 0.81])
g.set(xlim=(-0.41, 0.81))
g.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/expression_correlation_ridge_updated.png', dpi=500, bbox_inches='tight')
g.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/expression_correlation_ridge_updated.pdf', dpi=500, bbox_inches='tight')
plt.clf()

#Per severity
import seaborn as sns
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)}, font_scale=0.8)
# Define the data and variables
df = pd.DataFrame(pd.concat([wb.corr_mild, wb.corr_sev], axis=0))
df.columns = ["r"]
df["severity"] = ["mild"]*int(df.shape[0]/2)+["sev"]*int(df.shape[0]/2)
x_var = "r"
y_var = "severity"
means = df.groupby(y_var)[x_var].mean()
# Create the ridge plot
plt.figure(figsize=(12,6))
g = sns.FacetGrid(df, row=y_var, hue=y_var, aspect=10, height=0.5, palette=dict(zip(["mild", "sev"], ["royalblue", "blueviolet"])), row_order=["mild", "sev"])
g.map(sns.kdeplot, x_var, clip_on=False, fill=True, alpha=1, lw=1.5, bw_method=.25, clip=[-0.41, 0.81])
g.map(sns.kdeplot, x_var, clip_on=False, color="w", lw=2, bw_method=.25, clip=[-0.41, 0.81])
g.map(plt.axhline, y=0, lw=2, clip_on=False)
g.map(lambda **kwargs: plt.axvline(0, linestyle='--', color='tab:gray'), alpha=0.5, lw=1)
g.axes[0, 0].axvline(x=means["mild"], linestyle='--', color='black')
g.axes[1, 0].axvline(x=means["sev"], linestyle='--', color='black')
g.axes[0, 0].text(x=-0.397, y=2, s="COVID-19 Severity:\nAsymptomatic / Mild", color="royalblue", fontsize=8)
g.axes[1, 0].text(x=-0.397, y=1.7, s="COVID-19 Severity:\nMost Severe", color="blueviolet", fontsize=8)
g.set(yticks=[], ylabel=None, xlabel="Corr(mRNA, protein) in Pearson R")
g.set_titles("")
g.axes[1, 0].set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
g.axes[1, 0].set_ylabel("        Density", zorder=10, fontsize=10)
plt.xlim([-0.41, 0.81])
g.set(xlim=(-0.41, 0.81))
g.savefig('/Users/qingbowang/Desktop/taskforce_n1102/plots/expression_corr_persev_ridge.png', dpi=500, bbox_inches='tight')
plt.clf()

#stratify by %expression class
plt.rcParams.update({'font.size': 12})
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)}, font_scale=1)
fig, ax = plt.subplots(1, 5, sharex=True, sharey=True, figsize=(8, 5))
palette = sns.cubehelix_palette(5, start=2)
#titles = ["<1", "[1,10)","[10,25)","[25,50)","50$\leq$"]
titles = ["<1", "[1,5)","[5,25)","[25,50)","50$\leq$"]
i = 0
for b in wb.bin_detail.sort_values().unique():
    wbs = wb[wb.bin_detail==b]
    for row in wbs.index:
        y1, y2 = wbs.corr_mild[row], wbs.corr_sev[row]
        if y1>y2: c = "tab:blue"
        else: c = "tab:red"
        ax[i].plot([0, 1], [y1, y2], color=c, linestyle='-', linewidth=1, alpha=0.6, zorder=-2)
    ax[i].scatter([0]*len(wbs.corr_mild), wbs.corr_mild, color=palette[i], zorder=1, alpha=0.6)
    ax[i].scatter([0] * len(wbs.corr_mild), wbs.corr_mild, edgecolor="black", color="None", zorder=1, alpha=0.4) #for the edge
    ax[i].scatter([1]*len(wbs.corr_sev), wbs.corr_sev, color=palette[i], zorder=1, alpha=0.6)
    ax[i].scatter([1] * len(wbs.corr_sev), wbs.corr_sev, edgecolor="black", color="None", zorder=1, alpha=0.4)  # for the edge
    ax[i].set_xlabel("%Whole blood\nexpression: {0}".format(titles[i]), color=palette[i])
    ax[i].set_xticks([0,1],["Mild /\nAsymptomatic", "Most\nSevere"], rotation=60)
    i += 1
ax[0].set_ylabel("Corr(mRNA, protein) in Pearson R")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/ep_corr_per_severity_per_expression.png", dpi=500)#Overall
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/ep_corr_per_severity_per_expression.pdf", dpi=500)
plt.clf()


#Also plot only the essential stats
#i.e. mean pearson R for mild vs severe, overall and stratified by the %WB
plt.rcdefaults()
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
sns.reset_orig()
#sns.set(style="white", font_scale=1)
fig,ax=plt.subplots(figsize=(4.5,3.5))
sbd = pd.concat([wb[["bin_detail", "corr_mild"]].rename(columns={"corr_mild":"Corr(mRNA, protein)\nPearson R"}), wb[["bin_detail", "corr_sev"]].rename(columns={"corr_sev":"Corr(mRNA, protein)\nPearson R"})], axis=0)
sbd["Severity"] = ["Mild / Asymptomatic"]*wb.shape[0] + ["Most Severe"]*wb.shape[0]
sbd.rename(columns={"bin_detail":"%Whole Blood expression in GTEx"}, inplace=True)
sns.violinplot(data=sbd, x="%Whole Blood expression in GTEx", y="Corr(mRNA, protein)\nPearson R", hue="Severity",
               split=True, inner="quart", fill=False,
               palette={"Mild / Asymptomatic": "royalblue", "Most Severe": "blueviolet"})
plt.ylim([-0.41, .81])
ax.legend(loc='upper left', facecolor='white', bbox_to_anchor=(0.25,1.3))
#ax.set_xticklabels(["<1", "[1,10)","[10,25)","[25,50)","50$\leq$"])
ax.set_xticklabels(["<1", "[1,5)","[5,25)","[25,50)","50$\leq$"])
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/ep_corr_per_severity_per_expression_violin.png", dpi=500)#Overall
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/ep_corr_per_severity_per_expression_violin.pdf", dpi=500)#Overall
plt.clf()

#Scatter plots for individual gene, IL1RL1 = ENSG00000115602.16
ilm = dfm.loc[:,"ENSG00000115602.16"]
ilp = dfp.loc[:,"ENSG00000115602.16"]
meta = pd.read_csv("~/Desktop/taskforce_n1102/n1300/eqtl_n998.combined_covariates_idadded.txt", sep="\t", index_col=0)
meta = meta.iloc[-3:,:].T
palette = ["tab:gray", "lightsteelblue", "royalblue", "blueviolet"]
titles = ["Asymptomatic", "Mild", "Severe", "Most severe"]
fig, ax = plt.subplots(1, 4, sharex=True, sharey=True, figsize=(10,2.5))
for sev in [1,2,3,4]:
    x = ilm[meta.severity == sev]
    y = ilp[meta.severity == sev]
    r, p = stats.pearsonr(x, y)
    ax[sev-1].scatter(x, y, color=palette[sev-1], alpha=0.75, edgecolor=palette[sev-1])
    ax[sev-1].set_title(titles[sev-1], color=palette[sev-1])
    ax[sev-1].text(-3.2,2.75,"r={0}".format(r.round(4)), color=palette[sev-1], fontsize=14 )
    ax[sev-1].set_aspect('equal')
    ax[0].set_ylabel("Protein expression")
    ax[sev-1].set_xlabel("mRNA expression")
    ax[sev-1].axvline(x=0, color="black", linestyle="--", linewidth=0.5, zorder=-2)
    ax[sev-1].axhline(y=0, color="black", linestyle="--", linewidth=0.5, zorder=-2)
ax[0].set_xlim([-3.5,3.5])
ax[0].set_ylim([-3.5,3.5])
ax[0].set_xticks([-3,-2,-1,0,1,2,3])
ax[0].set_yticks([-3,-2,-1,0,1,2,3])
fig.suptitle('IL1RL1 expression per COVID-19 severity')
#plt.tight_layout()
plt.subplots_adjust(wspace=0.05, bottom=0.2, top=0.8)
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/IL1RL1_scatter.pdf", dpi=500)#Overall
plt.savefig("/Users/qingbowang/Desktop/taskforce_n1102/plots/IL1RL1_scatter.png", dpi=500)#Overall
plt.clf()
