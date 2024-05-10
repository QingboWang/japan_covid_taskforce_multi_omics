library(susieR)
library(data.table) #need to download later when the internet access is there
args = commandArgs(trailingOnly=TRUE)
ld = args[1]
z = args[2]
cs_out = args[3]
pip_out = args[4]
alpha_out = args[5]
n = as.integer(args[6])

R = as.matrix(fread(ld)) #check whether the format is still good with this fread.
print ("done reading ld in R")
print (Sys.time())
st = read.table(z, header=T)
st$z = st$beta / st$se

z_scores = st$z
print ("Starting SuSiE")
print (Sys.time())
fitted_rss = susie_rss(z_scores, R, n=n, L=5, refine=T, estimate_residual_variance = TRUE, max_iter = 1000) #yes do this refinement. Taking time is fine. And we set L=5 this time.
print ("Done SuSiE, starting writing results")
print (Sys.time())
#write the results
write.table(summary(fitted_rss)$cs , cs_out, quote=F, sep='\t')
st$pip = fitted_rss$pip
write.table(st[c("rsid","pip")], pip_out, quote=F, sep='\t')

write.table(fitted_rss$alpha, alpha_out, quote=F, sep='\t')

print ("Done writing")
print (Sys.time())

