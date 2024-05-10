#get the n-th gene name etc in the expression file

#print the num
echo "Starting $1-th line"

bedf=/Users/qingbowang/Desktop/taskforce_n1102/n1369.protein.expression.bed #says 1369 but I know it is 1384
l=$(head -n $1 $bedf | tail -n 1 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)

#prep. the vcf in cis-window
fn=/Users/qingbowang/Desktop/taskforce_n1102/n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_"$chr"_reIDed.vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/n1419_vcf_"$gn".tsv
#cut the file
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix -h $fn "$chr":$(($tss-1000000))-$(($tss+1000000)) > $fn_tmp
z_pqtl=/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n1384/"$gn"_p_fminput.z #z score file as a reference just to match the index
#create the LD matrix for n1384 pQTL:
now=$(date +"%T")
echo "Starting to create LD : $now"
python3 /Users/qingbowang/Desktop/taskforce_n1102/n1300/create_ld_mat_for_a_gene_pqtl.py $gn $fn_tmp $z_pqtl #Check that this is customized for n1384. yes it is now.
now=$(date +"%T")
echo "Done create LD : $now"
rm $fn_tmp #not needed anymore

#1. pQTL fine-mapping
now=$(date +"%T")
echo "Starting pQTL : $now"
#w/ FINEMAP
now=$(date +"%T")
echo "Starting to run FINEMAP : $now"
mkdir -p /Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/pqtl_n1384/"$gn"
cd /Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/pqtl_n1384/"$gn"
ln -s $z_pqtl ./"$gn"_p_fminput.z
ln -s /Users/qingbowang/Desktop/tmp/ld_covadj_pqtl_"$gn".ld ./"$gn"_p.ld
echo "z;ld;snp;config;cred;log;n_samples" > ./master
echo "${gn}_p_fminput.z;${gn}_p.ld;${gn}_p.snp;${gn}_p.config;${gn}_p.cred;${gn}_p.log;1384" >> ./master
~/Downloads/finemap_v1.3.1_MacOSX/finemap_v1.3.1_MacOSX n --sss --in-files ./master
now=$(date +"%T")
echo "Done FINEMAP $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/n1384_pqtl_finemap_reimputed_endlog.txt

#w/ SuSiE
now=$(date +"%T")
echo "Starting to run SuSiE : $now"
ld=/Users/qingbowang/Desktop/tmp/ld_covadj_pqtl_"$gn".ld
cs_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/pqtl_n1384/"$gn"/"$gn"_p_susie_cs.txt
pip_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/pqtl_n1384/"$gn"/"$gn"_p_susie_pip.tsv
alpha_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/pqtl_n1384/"$gn"/"$gn"_p_susie_alpha.tsv
n=1384
Rscript /Users/qingbowang/Desktop/codes/run_susie_for_a_gene.R $ld $z_pqtl $cs_out $pip_out $alpha_out $n
now=$(date +"%T")
echo "Done running SuSiE (l=$1) : $now"
echo "Done SuSiE $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/n1384_pqtl_susie_reimputed_endlog.txt

#remove the LD mat and subset vcf since it is so heavy
rm /Users/qingbowang/Desktop/tmp/ld_covadj_pqtl_"$gn".ld



