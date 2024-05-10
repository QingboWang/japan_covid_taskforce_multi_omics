#Collect the fm results
cd /Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary
gunzip -c /Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.mhcrem.bed.gz | cut -f4 > genes_ids.tsv
gene_names=genes_ids.tsv
gn_eg=ENSG00000172878.14 #just random gene name to begin with
#writing the header:
#eQTL
head -n 1 ../eqtl_n1019/"$gn_eg"/"$gn_eg"_e.snp | cut -d " " -f2,11 > n1019_eqtl_finemap_pip0001.txt
echo "gene_id\t n_var_below_thres" > n1019_eqtl_finemap_n_below_pip0001.txt
head -n 1 ../eqtl_n1019/"$gn_eg"/"$gn_eg"_e_susie_pip.tsv > n1019_eqtl_susie_pip0001.txt
echo "gene_id\t n_var_below_thres" > n1019_eqtl_susie_n_below_pip0001.txt

#count
i=0
while IFS= read -r line
do
    #FINEMAP
    tail -n +2 ../eqtl_n1019/"$line"/"$line"_e.snp | cut -d " " -f2,11 | awk -F " " '$2 >= 0.001' >> n1019_eqtl_finemap_pip0001.txt
    n_var_belowthres=`cat ../eqtl_n1019/"$line"/"$line"_e.snp | awk -F " " '$11 < 0.001' | wc -l`
    echo "$line\t$n_var_belowthres" >> n1019_eqtl_finemap_n_below_pip0001.txt
    #SuSiE:
    tail -n +2 ../eqtl_n1019/"$line"/"$line"_e_susie_pip.tsv | cut -f2,3 | awk '$2 >= 0.001' >> n1019_eqtl_susie_pip0001.txt
    n_var_belowthres=`cat ../eqtl_n1019/"$line"/"$line"_e_susie_pip.tsv | awk '$3 < 0.001' | wc -l`
    echo "$line\t$n_var_belowthres" >> n1019_eqtl_susie_n_below_pip0001.txt
    i=$(($i+1))
    if [ $(($i%100)) -eq 0 ]
    then
      now=$(date +"%T")
      echo "done $i $now"
    fi
done < "$gene_names"

cd /Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/summary
cut -f4 /Users/qingbowang/Desktop/taskforce_n1102/n1369.protein.expression.bed > pqtl_genes_ids.tsv #says 1369 but I know it is 1384
gene_names=pqtl_genes_ids.tsv
gn_eg=ENSG00000188157.15 #just random gene name to begin with
#writing the header:
head -n 1 ../pqtl_n1384/"$gn_eg"/"$gn_eg"_p.snp | cut -d " " -f2,11 > n1384_pqtl_finemap_pip0001.txt
echo "gene_id\t n_var_below_thres" > n1384_pqtl_finemap_n_below_pip0001.txt
head -n 1 ../pqtl_n1384/"$gn_eg"/"$gn_eg"_p_susie_pip.tsv > n1384_pqtl_susie_pip0001.txt
echo "gene_id\t n_var_below_thres" > n1384_pqtl_susie_n_below_pip0001.txt

gene_names=pqtl_genes_ids.tsv
while IFS= read -r line
do
  mv ../pqtl_n1384/"$line"/"$line"_e.snp ../pqtl_n1384/"$line"/"$line"_p.snp
done < "$gene_names"

#count
i=0
while IFS= read -r line
do
    #FINEMAP
    tail -n +2 ../pqtl_n1384/"$line"/"$line"_p.snp | cut -d " " -f2,11 | awk -F " " '$2 >= 0.001' >> n1384_pqtl_finemap_pip0001.txt
    n_var_belowthres=`cat ../pqtl_n1384/"$line"/"$line"_p.snp | awk -F " " '$11 < 0.001' | wc -l`
    echo "$line\t$n_var_belowthres" >> n1384_pqtl_finemap_n_below_pip0001.txt
    #SuSiE:
    tail -n +2 ../pqtl_n1384/"$line"/"$line"_p_susie_pip.tsv | cut -f2,3 | awk '$2 >= 0.001' >> n1384_pqtl_susie_pip0001.txt
    n_var_belowthres=`cat ../pqtl_n1384/"$line"/"$line"_p_susie_pip.tsv | awk '$3 < 0.001' | wc -l`
    echo "$line\t$n_var_belowthres" >> n1384_pqtl_susie_n_below_pip0001.txt
    i=$(($i+1))
    if [ $(($i%100)) -eq 0 ]
    then
      now=$(date +"%T")
      echo "done $i $now"
    fi
done < "$gene_names"