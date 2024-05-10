#RNA QC metrics from RNAseQC:
import pandas as pd
qc1 = pd.read_csv("~/Downloads/rnaseqc_metrics_n500.tsv", sep='\t', index_col=0)
qc2 = pd.read_csv("~/Downloads/rnaseqc_metrics_n602.tsv", sep='\t', index_col=0)
qc1["batch"] = "early"
qc2["batch"] = "late"
qc = pd.concat([qc1, qc2.iloc[:,1:]])

#take a look at different metrics
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
mu = qc.Mapped.mean()
sigma = qc.Mapped.std()
plt.figure(figsize=(4.5,3))
p = qc.Mapped.hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    if ((rectangle.get_x() < 0.5*10**8) | (rectangle.get_x() > 3*10**8)):
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.5*10**8, linestyle="--", linewidth=1, color="black")
plt.axvline(x=3*10**8, linestyle="--", linewidth=1, color="black")
plt.axvline(x=mu, linestyle=":", linewidth=1, color="black")
plt.xlabel("Number of mapped reads")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_mappedreads_n1102_customthres.png", dpi=500)
plt.clf()

mu = qc["Mapping Rate"].mean()
sigma = qc["Mapping Rate"].std()
plt.figure(figsize=(4.5,3))
p = qc["Mapping Rate"].hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    #if ( (rectangle.get_x() < mu-2*sigma)|(rectangle.get_x() > mu+2*sigma) ):
    if ((rectangle.get_x() < 0.97)):#high -> OK
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.97, linestyle="--", linewidth=1, color="black")
plt.axvline(x=mu, linestyle=":", linewidth=1, color="black")
plt.xlabel("Mapping rate")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_maprate_n1102_custom.png", dpi=500)
plt.clf()

mu = qc["Intergenic Rate"].mean()
sigma = qc["Intergenic Rate"].std()
plt.figure(figsize=(4.5,3))
p = qc['Intergenic Rate'].hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    #if ( (rectangle.get_x() < mu-2*sigma)|(rectangle.get_x() > mu+2*sigma) ):
    if ((rectangle.get_x() > 0.05)): #Low -> OK
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.05, linestyle="--", linewidth=1, color="black")
plt.axvline(x=mu, linestyle=":", linewidth=1, color="black")
plt.xlabel("Intergenic rate")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_intgenic_n1102_custom.png", dpi=500)
plt.clf()

mu = qc["Base Mismatch Rate"].mean()
sigma = qc["Base Mismatch Rate"].std()
plt.figure(figsize=(4.5,3))
p = qc['Base Mismatch Rate'].hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    #if ( (rectangle.get_x() < mu-2*sigma)|(rectangle.get_x() > mu+2*sigma) ):
    if ((rectangle.get_x() > 0.005)): #Low -> OK
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.005, linestyle="--", linewidth=1, color="black")
plt.axvline(x=mu, linestyle=":", linewidth=1, color="black")
plt.xlabel("Base mismatch rate")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_basemismatch_n1102_custom.png", dpi=500)
plt.clf()

mu = qc["rRNA rate"].mean()
sigma = qc["rRNA rate"].std()
plt.figure(figsize=(4.5,3))
p = qc['rRNA rate'].hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    #if ( (rectangle.get_x() < mu-2*sigma)|(rectangle.get_x() > mu+2*sigma) ):
    if ((rectangle.get_x() > 0.05)): #Low -> OK
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=0.05, linestyle="--", linewidth=1, color="black")
plt.axvline(x=mu, linestyle=":", linewidth=1, color="black")
plt.xlabel("Ribosomal RNA rate")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_rrna_n1102_custom.png", dpi=500)
plt.clf()

D_i = pd.read_csv("~/Downloads/cnt_n1102_corr_D_i.tsv", sep='\t', index_col=0, squeeze=True)
mu = D_i.mean()
sigma = D_i.std()
plt.figure(figsize=(4.5,3))
p = D_i.hist(bins=50, grid=False, color="tab:blue")
for rectangle in p.patches:
    #if ( (rectangle.get_x() < mu-2*sigma)|(rectangle.get_x() > mu+2*sigma) ):
    if ((rectangle.get_x() < -15)): #High -> OK
        rectangle.set_facecolor('tab:gray')
plt.axvline(x=-15, linestyle="--", linewidth=1, color="black") #これも2sdにしてみるか
plt.xlabel("Intersample correlation deviation ($D_i$)")
plt.ylabel("Number of samples")
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/qc_di_n1102_custom.png", dpi=500)
plt.clf()

#check each filtering criteria:
filtered_by_n_mappedread = ((qc.Mapped<0.5*10**8)|(qc.Mapped>3*10**8))
filtered_by_mapping_rate = ((qc["Mapping Rate"]<0.97))
filtered_by_intergenic_rate = (qc["Intergenic Rate"]>0.05)
filtered_by_base_mismatch_rate = (qc["Base Mismatch Rate"]>0.005)
filtered_by_rrna_rate = (qc["rRNA rate"]>0.05)
filtered_by_intersample_corr = (D_i<-15)
filtered_by_intersample_corr.name = "Intersample_corr"
filt = pd.concat([filtered_by_n_mappedread, filtered_by_mapping_rate,
           filtered_by_intergenic_rate, filtered_by_base_mismatch_rate,
           filtered_by_rrna_rate, filtered_by_intersample_corr], axis=1)
filt.value_counts().sort_index()
filt.value_counts().sort_index().to_csv("~/Downloads/rnaseq_n1102_filter_custom.tsv", sep="\t")
