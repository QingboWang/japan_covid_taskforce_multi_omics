##a function to calc. in-sample cov-adjusted LD mat for +-1mb of each of the gene tss

""" #This is e.g. how to run this python code
l=$(gunzip -c /Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.bed.gz | head -n $n | tail -n 1 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)
fn=/Users/qingbowang/Desktop/taskforce_n1102//n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_"$chr"_reIDed.vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/n1419_vcf_"$gn".tsv
#cut the file
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix -h $fn "$chr":$(($tss-1000000))-$(($tss+1000000)) > $fn_tmp
z_eqtl=/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/eqtl/"$gn"_e_fminput.z #z score file as a reference just to match the index
#z_pqtl=/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl/"$gn"_p_fminput.z #z score file as a reference just to match the index
#z_esqtl=/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/esqtl/"$gn"_es_fminput.z #z score file as a reference just to match the index
#z_psqtl=/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/psqtl/"$gn"_ps_fminput.z #z score file as a reference just to match the index
python3 /Users/qingbowang/Desktop/taskforce_n1102/n1300/create_ld_mat_for_a_gene.py $gn $fn_tmp $z_eqtl
rm $fn_tmp
"""

#input = the file, and also the gene z file as a reference
import sys
import os
from os.path import exists
import io
import numpy as np
import pandas as pd
import time as tm
gene = sys.argv[1]
df_dir = sys.argv[2]
z_eqtl_dir = sys.argv[3]
#z_pqtl_dir = sys.argv[3] #If pQTL use this

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

df = read_vcf(df_dir)
df.index = df.ID #since it is already well re-IDed
#set the order same (which should be, by default)
z = pd.read_csv(z_eqtl_dir, sep=' ')
#z = pd.read_csv(z_pqtl_dir, sep=' ')
df = df.loc[z.rsid.str.split("_ENSG").str[0],:]

#and turn it to dosage matrix;
X0 = df.iloc[:,9:].applymap(lambda x: float(x.split(":")[1])).T
print ("done getting dosage mat for gene{0}: {1}".format(gene, tm.ctime()))

#and filter to n=1019 samples:
s = pd.read_csv("~/Desktop/taskforce_n1102/taskforce_n1019_sex_age_sev.tsv", sep='\t', index_col=0).T
#to n=1384 samples for the full pQTL call
#s = pd.read_csv("~/Desktop/taskforce_n1102/n1300/taskforce_n1384_forprotein_sex_age_sev.tsv", sep='\t', index_col=0).T
#X0 = X0.loc[:, s.index]
X0 = X0.loc[s.index, :] #since we transposed

#print ("First is for eQTL, starting creating X^tX/N, {0}".format(tm.ctime()))
print ("for pQTL, starting creating X^tX/N, {0}".format(tm.ctime()))
X = X0.apply(lambda l: (l - l.mean())/l.std(), axis=0)
M = pd.read_csv("~/Desktop/taskforce_n1102/n1019_eqtl_M_for_covadj.tsv", sep="\t", index_col=0) #TBD
#M = pd.read_csv("~/Desktop/taskforce_n1102/n1300/n1384_pqtl_M_for_covadj.tsv", sep="\t", index_col=0) #TBD
#get the adj matrix
M = M.loc[X.index, X.index] #make sure the order is the same
X_adj = M @ X
n = X_adj.shape[0]#number of samples
R = X_adj.T @ X_adj / n
print ("done calculating LD mat for gene{0}: {1}".format(gene, tm.ctime()))
print (R.shape)
print (R.iloc[:5,:5])
#and write as .ld:
R = np.array(R)
print ("starting writing, {0}".format(tm.ctime()))
np.savetxt("/Users/qingbowang/Desktop/tmp/ld_covadj_eqtl_{0}.ld".format(gene), R, delimiter=' ', fmt='%4f')
#np.savetxt("/Users/qingbowang/Desktop/tmp/ld_covadj_pqtl_{0}.ld".format(gene), R, delimiter=' ', fmt='%4f')
print ("done writing, {0}".format(tm.ctime()))