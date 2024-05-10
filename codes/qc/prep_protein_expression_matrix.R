library(OlinkAnalyze)
library(openxlsx)
library(dplyr)

setwd("~/Desktop/taskforce_n1102/olink_data/")

#read NPX data, three  batches, already filtered to CT samples:
df_Q00984 = read_NPX("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-00984_NPX_warningfilt_filtered_to_CT.csv")
df_Q01085 = read_NPX("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-01085_NPX_warningfilt_filtered_to_CT.csv")
df_Q01086 = read_NPX("~/Desktop/taskforce_n1102/olink_data/cleaned/Q-01086_NPX_warningfilt_filtered_to_CT.csv")

#Subset normalization, adding 01085:
normed = olink_normalization(df1 = df_Q00984[,colnames(df_Q01085)], #since this one has more colnames
                        df2 = df_Q01085,
                        overlapping_samples_df1 = unique(df_Q00984$SampleID),
                        overlapping_samples_df2 = unique(df_Q01085$SampleID))
#Bridging normalization using 16 overlapping samples, adding 01086:
intersect(normed[,"Variable.2..Subject."], df_Q01086[,"Variable.2..Subject."]) #16 samples that intersects
ct_ids_int = intersect(normed$SampleID, df_Q01086$SampleID) #is the SampleID in olink
#and use that to perform bridging normalization
normed2 = olink_normalization(df1 = normed[,colnames(df_Q01086)],
                        df2 = df_Q01086,
                        overlapping_samples_df1 = ct_ids_int)
normed2 = normed2 %>% select(-Project) %>% select(-Adj_factor) #removing redundant columns
write.table(normed2, "~/Desktop/taskforce_n1102/olink_data/n1300_normed_npx_20230114.tsv", sep='\t')


