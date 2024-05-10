import pandas as pd

eq = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.bed.gz", sep="\t")
pq = pd.read_csv("~/Desktop/taskforce_n1102/n1369.protein.expression.bed", sep='\t')
#intersect the genes
eq.index = eq.gene_id
pq.index = pq.gene_id
intsct = eq.index.intersection(pq.index)
eq.loc[intsct,:]
pq.loc[intsct,:]
sum(eq.loc[intsct,"start"]!=pq.loc[intsct,"start"]) #is zero. Perfect.
#make sure that the .version in ENSG is not the problem - done.
print (len(intsct), len(eq.index.str.split("\\.").str[0].intersection(pq.index.str.split("\\.").str[0])))
#intersect the individuals:
intsct_ct = eq.columns[4:].intersection(pq.columns[4:])
eq.loc[intsct,intsct_ct] #998 samples in the intersection

#removing mhc genes
left = 25726063
right = 33400644
df_mhc = eq[(eq["#chr"]=="chr6")&(eq.start>left-10**6)&(eq.start<right+10**6)]
mhc_genes = df_mhc.index.get_level_values(0)
intsct = intsct.difference(mhc_genes)

#Now filter to intersections:
eqfilt = eq.loc[intsct,intsct_ct]
pqfilt = pq.loc[intsct,intsct_ct]
#re-perform inverse normal transformation within these samples:
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
eqfilt_invnorm = inverse_normal_transform(eqfilt)
pqfilt_invnorm = inverse_normal_transform(pqfilt)
#Save
eqfilt.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/mrna_expression_intersection.tsv.gz", sep="\t", compression="gzip")
pqfilt.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/protein_expression_intersection.tsv.gz", sep="\t", compression="gzip")
eqfilt_invnorm.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/mrna_expression_intersection_invnormed.tsv.gz", sep="\t", compression="gzip")
pqfilt_invnorm.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/protein_expression_intersection_invnormed.tsv.gz", sep="\t", compression="gzip")

#Also re-annotate the index and output for eQTL and pQTL call in n=998 samples
eq = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/mrna_expression_intersection_invnormed.tsv.gz", sep="\t", compression="gzip", index_col=0)
pq = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/protein_expression_intersection_invnormed.tsv.gz", sep="\t", compression="gzip", index_col=0)
pqorig = pd.read_csv("~/Desktop/taskforce_n1102/n1369.protein.expression.bed", sep='\t')
pqorig.index = pqorig.gene_id
eq = eq.join(pqorig.iloc[:,:3], how="left").set_index(list(pqorig.columns[:3]), append=True)
pq = pq.join(pqorig.iloc[:,:3], how="left").set_index(list(pqorig.columns[:3]), append=True)
eq.reorder_levels([1,2,3,0]).sort_index().to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/mrna_expression_intersection_invnormed.bed", sep="\t")
pq.reorder_levels([1,2,3,0]).sort_index().to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/protein_expression_intersection_invnormed.bed", sep="\t")


#Now do the linear regression to regress out the other:
eq = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/mrna_expression_intersection_invnormed.tsv.gz", sep="\t", compression="gzip", index_col=0)
pq = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/protein_expression_intersection_invnormed.tsv.gz", sep="\t", compression="gzip", index_col=0)

#to do so, first make sure that a linear relationship is plausible (e.g. is log-linear better?)
#i.e. calculate the pearson correlation, and inspect the top and bottom 10 just as an example:
from scipy import stats
corr = pd.Series(eq.index).apply(lambda x: stats.pearsonr(eq.loc[x,:], pq.loc[x,:])[0])
corr.index = eq.index
corr.sort_values(inplace=True, ascending=False)
#plot top 10
for gene in corr.index[:10]:
    plt.scatter(eq.loc[gene,:], pq.loc[gene,:])
    plt.title(gene)
    plt.xlabel("mRNA")
    plt.ylabel("Protein")
    plt.show() #great, linear makes sense, although dropout seems to be a problem...
#e.g. worst 10
for gene in corr.index[-10:]:
    plt.scatter(eq.loc[gene,:], pq.loc[gene,:])
    plt.title(gene)
    plt.xlabel("mRNA")
    plt.ylabel("Protein")
    plt.show() #Also linear
#e.g. zero correlation:
corrdiffzero = (abs(corr)).sort_values(ascending=True)
for gene in corrdiffzero.index[:10]:
    plt.scatter(eq.loc[gene,:], pq.loc[gene,:])
    plt.title(gene)
    plt.xlabel("mRNA")
    plt.ylabel("Protein")
    plt.show() #Literally zero correlation yeah.
#conclusion: linear is great
def regress_out_vec(x,y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return (y - x*slope)
def regress_out_pergene(X, Y):#row = gene, columns = samples. returns Y where X is regressed out
    Y_adj = pd.Series(np.arange(Y.shape[0])).apply(lambda i: regress_out_vec(X.iloc[i,:], Y.iloc[i,:]))
    Y_adj.index = Y.index
    return (Y_adj)
rna = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/mrna_expression_intersection_invnormed.tsv.gz", sep="\t", compression="gzip", index_col=0)
pro = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/protein_expression_intersection_invnormed.tsv.gz", sep="\t", compression="gzip", index_col=0)
rna_adjusted = regress_out_pergene(pro, rna)
protein_adjusted = regress_out_pergene(rna, pro)
#re-invnorm:
es_invnorm = inverse_normal_transform(rna_adjusted)
ps_invnorm = inverse_normal_transform(protein_adjusted)
#add the gene info and save as bed to plug into es and psQTL call
pqorig = pd.read_csv("~/Desktop/taskforce_n1102/n1369.protein.expression.bed", sep='\t')
pqorig.index = pqorig.gene_id
es_invnorm = es_invnorm.join(pqorig.iloc[:,:3], how="left").set_index(list(pqorig.columns[:3]), append=True)
ps_invnorm = ps_invnorm.join(pqorig.iloc[:,:3], how="left").set_index(list(pqorig.columns[:3]), append=True)
#Save
es_invnorm.reorder_levels([1,2,3,0]).sort_index().to_csv("~/Desktop/taskforce_n1102/n1300/n998_mrna.protein_regressedout.bed", sep='\t')
ps_invnorm.reorder_levels([1,2,3,0]).sort_index().to_csv("~/Desktop/taskforce_n1102/n1300/n998_protein.mrna_regressedout.bed", sep='\t')




