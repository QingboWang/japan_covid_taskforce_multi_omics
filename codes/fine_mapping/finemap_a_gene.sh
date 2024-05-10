#get the n-th gene name etc in the expression file

#print the num
echo "Starting $1-th line"

l=$(gunzip -c /Users/qingbowang/Desktop/taskforce_n1102/n1019.expression.mhcrem.bed.gz | head -n $1 | tail -n 1 | cut -f1-4)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f3)
gn=$(echo $l | cut -f4)

#prep. the vcf in cis-window
fn=/Users/qingbowang/Desktop/taskforce_n1102/n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_"$chr"_reIDed.vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/n1419_vcf_"$gn".tsv
#cut the file
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix -h $fn "$chr":$(($tss-1000000))-$(($tss+1000000)) > $fn_tmp
z_eqtl=/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/eqtl_n1019/"$gn"_e_fminput.z #z score file as a reference just to match the index
#create the LD matrix for n1019 eQTL:
now=$(date +"%T")
echo "Starting to create LD : $now"
python3 /Users/qingbowang/Desktop/taskforce_n1102/n1300/create_ld_mat_for_a_gene.py $gn $fn_tmp $z_eqtl
now=$(date +"%T")
echo "Done create LD : $now"
rm $fn_tmp #not needed anymore

#1. eQTL fine-mapping
now=$(date +"%T")
echo "Starting eQTL : $now"
#w/ FINEMAP
now=$(date +"%T")
echo "Starting to run FINEMAP : $now"
mkdir -p /Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/eqtl_n1019/"$gn"
cd /Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/eqtl_n1019/"$gn"
ln -s $z_eqtl ./"$gn"_e_fminput.z
ln -s /Users/qingbowang/Desktop/tmp/ld_covadj_eqtl_"$gn".ld ./"$gn"_e.ld
echo "z;ld;snp;config;cred;log;n_samples" > ./master
echo "${gn}_e_fminput.z;${gn}_e.ld;${gn}_e.snp;${gn}_e.config;${gn}_e.cred;${gn}_e.log;1019" >> ./master
~/Downloads/finemap_v1.3.1_MacOSX/finemap_v1.3.1_MacOSX n --sss --in-files ./master
now=$(date +"%T")
echo "Done FINEMAP $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/n1019_eqtl_finemap_reimputed_endlog.txt

#w/ SuSiE
now=$(date +"%T")
echo "Starting to run SuSiE : $now"
ld=/Users/qingbowang/Desktop/tmp/ld_covadj_eqtl_"$gn".ld
cs_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/eqtl_n1019/"$gn"/"$gn"_e_susie_cs.txt
pip_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/eqtl_n1019/"$gn"/"$gn"_e_susie_pip.tsv
alpha_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/eqtl_n1019/"$gn"/"$gn"_e_susie_alpha.tsv
n=1019
Rscript /Users/qingbowang/Desktop/codes/run_susie_for_a_gene.R $ld $z_eqtl $cs_out $pip_out $alpha_out $n
now=$(date +"%T")
echo "Done running SuSiE (l=$1) : $now"
echo "Done SuSiE $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/n1019_eqtl_susie_reimputed_endlog.txt

#remove the LD mat and subset vcf since it is so heavy
rm /Users/qingbowang/Desktop/tmp/ld_covadj_eqtl_"$gn".ld



