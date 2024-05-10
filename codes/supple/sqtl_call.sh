#Run fastqtl, per chr:
cd /home/qwang/n1019_splice
prefix=leafcutter_n1019
for i in {1..22} X
do
vcf=../n1419_vcf_hg38/taskforce_n1419_imputed_hg38_sorted_chr${i}_reIDed.vcf.gz
outname=leafcutter_n1019_redo.chr${i}
echo "start $i"
singularity exec -e ../gtex_eqtl-V8.simg \
/bin/bash -c "/opt/fastqtl/python/run_FastQTL_threaded.py ${vcf} leafcutter_n1019_perind.counts.gz.qqnorm_chr${i}.cleaned.bed.gz ${outname} \
    --covariates leafcutter_n1019_redo.combined_covariates_idadded.txt \
    --window 1e5 --chunks 8 --threads 8 --ma_sample_threshold 1"
done

#Filter by min(p) and mac:
for i in {1..22}; do
    python ~/Desktop/codes/filter_intron_mac2.py $i &
done
#Where the python code is as below:
```filter_intron_mac2.py
import pandas as pd
import numpy as np
import sys
import time as tm
chr = sys.argv[1]
s = pd.read_csv("~/Desktop/taskforce_n1102/n1300/sqtl_sumstats/leafcutter_n1019_redo.chr{0}.allpairs.txt.gz".format(chr), sep='\t')
print ("chr{0}".format(chr))
print ("n bef")
print(s.shape[0])
s = s[s.ma_count>2]
print ("n(vars) after: ")
print(s.shape[0])
s.to_csv("~/Desktop/taskforce_n1102/n1300/sqtl_sumstats/sqtl_n1019_redo.chr{0}.allpairs.mac2.txt.gz".format(chr),sep='\t', compression="gzip", index=False)
print ("done chr{0}, {1}".format(chr, tm.ctime()))
#Then filter by pval
print ("Then filter by pval")
print ("n bef")
print(s.shape[0])
intron_level_min_pval = s.groupby("gene_id").pval_nominal.min()
s.index = s.gene_id
s = s.join(intron_level_min_pval, rsuffix="_minperintron")
s = s[s.pval_nominal_minperintron<5e-8]
del s['pval_nominal_minperintron']
print ("n(vars) after: ")
print(s.shape[0])
s.to_csv("~/Desktop/taskforce_n1102/n1300/sqtl_sumstats/sqtl_n1019_redo.chr{0}.allpairs.mac2.intron5em8.txt.gz".format(chr),sep='\t', compression="gzip", index=False)
intron_level_min_pval.to_csv("~/Desktop/taskforce_n1102/n1300/sqtl_sumstats/intron_level_min_pval_chr{0}.tsv".format(chr), sep="\t")
print ("done chr{0}, {1}".format(chr, tm.ctime()))
```

#prepare the z files for fine-mapping:
for i in {1..22}; do
    python ~/Desktop/codes/write_sqtl_z.py $i &
done
#Where the python code is as below:
```write_sqtl_z.py
import time as tm
import gzip
import sys
output_path = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/sqtl/" #done: mkdir
chr = sys.argv[1]
i = 0
j = 0
input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/sqtl_sumstats/sqtl_n1019_redo.chr{0}.allpairs.mac2.intron5em8.txt.gz".format(chr)
outF = open("dummy_file.txt", "w") #dummy to close the file
with gzip.open(input_file) as f:
            header = " ".join(["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]) #fm input style
            gene_now = "" #the gene we are writing now
            for line in f:
                gene = line.decode("utf-8").split("\t")[0] #actually this is a cluster ID but fine..
                if gene != gene_now: #starting a new gene
                        outF.close()
                        gene_now = gene
                        print ("starting {0} on chr {1}, {2}th line, {3} th cluster, {4}".format(gene, chr, i, j, tm.ctime()))
                        outF = open("{0}{1}_s_fminput.z".format(output_path, gene), "w")
                        bytes = outF.write(header+"\n")
                        j += 1
                else:
                    outF = open("{0}{1}_s_fminput.z".format(output_path, gene), "a")
                    #parse for the output finemap style
                    parts = line.decode("utf-8").split("\t")
                    rsid = parts[1] + "_" + parts[0] #variant-intron
                    chromosome = parts[1].split(":")[0]
                    position = parts[1].split(":")[1]
                    allele1 = parts[1].split(":")[2]
                    allele2 = parts[1].split(":")[3].split("_")[0]
                    mac = parts[4]
                    maf = parts[5]
                    beta = parts[7]
                    se = parts[8]
                    towrite = (" ").join([rsid, chromosome, position, allele1, allele2, maf, beta, se])
                    #if float(mac)>1: #skip mac0 and 1 -> this is already done so no need
                    bytes = outF.write(towrite)
                i += 1
```
