import pandas as pd
import time as tm
#Just doing in-sample mac>2 filtering and write

for chr in list(range(1,22+1))+["X"]:
    p = pd.read_csv("~/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.txt.gz".format(chr), sep='\t')
    e = pd.read_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.txt.gz".format(chr), sep='\t')
    es = pd.read_csv("~/Desktop/taskforce_n1102/n1300/esqtl_sumstats/esqtl_n998.chr{0}.allpairs.txt.gz".format(chr), sep='\t')
    ps = pd.read_csv("~/Desktop/taskforce_n1102/n1300/psqtl_sumstats/psqtl_n998.chr{0}.allpairs.txt.gz".format(chr), sep='\t')
    print ("(e,p,es,ps)")
    print("({0},{1},{2},{3})".format(e.shape[0], p.shape[0], es.shape[0], ps.shape[0]))
    p = p[p.ma_count>2]
    e = e[e.ma_count>2] #Now the threshold is 2 not 1 (i.e. we need at least 3 minor allele count)
    ps = ps[ps.ma_count>2]
    es = es[es.ma_count>2]
    print ("n(vars) after: ")
    print ("(e,p,es,ps)")
    print("({0},{1},{2},{3})".format(e.shape[0], p.shape[0], es.shape[0], ps.shape[0]))
    p.to_csv("~/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t', compression="gzip", index=False)
    e.to_csv("/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr), sep='\t', compression="gzip", index=False)
    es.to_csv("~/Desktop/taskforce_n1102/n1300/esqtl_sumstats/esqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr),sep='\t', compression="gzip", index=False)
    ps.to_csv("~/Desktop/taskforce_n1102/n1300/psqtl_sumstats/psqtl_n998.chr{0}.allpairs.mac2.txt.gz".format(chr),sep='\t', compression="gzip", index=False)
    print ("done chr{0}, {1}".format(chr, tm.ctime()))







