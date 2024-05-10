#print the num
echo "Starting $1-th line"

l=$(head -n $1 /Users/qingbowang/Desktop/taskforce_n1102/n1300/sqtl_sumstats/introns_to_finemap_mhcrem.tsv | tail -n 1)
chr=$(echo $l | cut -f3)
tss=$(echo $l | cut -f4) #actually this is cluster center position, but fine
gn=$(echo $l | cut -f1) #actually this is intron cluster name, but fine

#prep. the vcf in the narrowed cis-window (+-0.1Mb)
fn=/Users/qingbowang/Desktop/taskforce_n1102/n1300/vcf_hg38/taskforce_n1419_imputed_hg38_sorted_"$chr"_reIDed.vcf.gz
fn_tmp=/Users/qingbowang/Desktop/tmp/n1419_vcf_"$gn".tsv
#cut the file
/Users/qingbowang/samtools-1.13/htslib-1.13/tabix -h $fn "$chr":$(($tss-100000))-$(($tss+100000)) > $fn_tmp #Just 0.1Mb around cluster center
z_sqtl=/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/sqtl/"$gn"_s_fminput.z #z score file as a reference just to match the index
#create the LD matrix:
now=$(date +"%T")
echo "Starting to create LD : $now"
python3 /Users/qingbowang/Desktop/taskforce_n1102/n1300/create_ld_mat_for_an_intron.py $gn $fn_tmp $z_sqtl
now=$(date +"%T")
echo "Done create LD : $now"
rm $fn_tmp #not needed anymore

#w/ FINEMAP
now=$(date +"%T")
echo "Starting to run FINEMAP : $now"
mkdir -p /Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/sqtl/"$gn"
cd /Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/sqtl/"$gn"
ln -s $z_sqtl ./"$gn"_s_fminput.z
ln -s /Users/qingbowang/Desktop/tmp/ld_covadj_sqtl_"$gn".ld ./"$gn"_s.ld
echo "z;ld;snp;config;cred;log;n_samples" > ./master
echo "${gn}_s_fminput.z;${gn}_s.ld;${gn}_s.snp;${gn}_s.config;${gn}_s.cred;${gn}_s.log;1019" >> ./master
~/Downloads/finemap_v1.3.1_MacOSX/finemap_v1.3.1_MacOSX n --sss --in-files ./master
now=$(date +"%T")
echo "Done FINEMAP $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/n1019_sqtl_redo_finemap_endlog.txt

#w/ SuSiE
now=$(date +"%T")
echo "Starting to run SuSiE : $now"
ld=/Users/qingbowang/Desktop/tmp/ld_covadj_sqtl_"$gn".ld
cs_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/sqtl/"$gn"/"$gn"_s_susie_cs.txt
pip_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/sqtl/"$gn"/"$gn"_s_susie_pip.tsv
alpha_out=/Users/qingbowang/Desktop/taskforce_n1102/n1300/fmoutputs/sqtl/"$gn"/"$gn"_s_susie_alpha.tsv
n=1019
Rscript /Users/qingbowang/Desktop/codes/run_susie_for_a_gene.R $ld $z_sqtl $cs_out $pip_out $alpha_out $n
now=$(date +"%T")
echo "Done running SuSiE (l=$1) : $now"
echo "Done SuSiE $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/n1019_sqtl_redo_susie_endlog.txt

#remove the LD mat and subset vcf since it is so heavy
rm /Users/qingbowang/Desktop/tmp/ld_covadj_sqtl_"$gn".ld