import time as tm
import gzip

#write eQTL stats, n1019, for all genes:
output_path = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/eqtl_n1019/" #done: mkdir
i = 0
j = 0
for chr in list(range(1,22+1))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/eqtl_sumstats/n1419.n1019.chr{0}.allpairs.mac2.txt.gz".format(chr)
    outF = open("dummy_file.txt", "w") #dummy to close the file
    with gzip.open(input_file) as f:
        header = " ".join(["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]) #fm input style
        gene_now = "" #the gene we are writing now
        for line in f:
            gene = line.decode("utf-8").split("\t")[0]
            if gene != gene_now: #starting a new gene
                    outF.close()
                    gene_now = gene
                    print ("starting {0}, {1}th line, {2} th gene, {3}".format(gene, i, j, tm.ctime()))
                    outF = open("{0}{1}_e_fminput.z".format(output_path, gene), "w")
                    bytes = outF.write(header+"\n")
                    j += 1
            else:
                outF = open("{0}{1}_e_fminput.z".format(output_path, gene), "a")
                #parse for the output finemap style
                parts = line.decode("utf-8").split("\t")
                rsid = parts[1] + "_" + parts[0] #variant-gene
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


#write pQTL stats, n1384, for all genes:
output_path = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/zfiles/pqtl_n1384/" #done: mkdir
i = 0
j = 0
for chr in list(range(1,22+1))+["X"]:
    input_file = "/Users/qingbowang/Desktop/taskforce_n1102/n1300/pqtl_sumstats/n1384.protein.chr{0}.allpairs.mac2.txt.gz".format(chr)
    outF = open("dummy_file.txt", "w") #dummy to close the file
    with gzip.open(input_file) as f:
        header = " ".join(["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]) #fm input style
        gene_now = "" #the gene we are writing now
        for line in f:
            gene = line.decode("utf-8").split("\t")[0]
            if gene != gene_now: #starting a new gene
                    outF.close()
                    gene_now = gene
                    print ("starting {0}, {1}th line, {2} th gene, {3}".format(gene, i, j, tm.ctime()))
                    outF = open("{0}{1}_p_fminput.z".format(output_path, gene), "w")
                    bytes = outF.write(header+"\n")
                    j += 1
            else:
                outF = open("{0}{1}_p_fminput.z".format(output_path, gene), "a")
                #parse for the output finemap style
                parts = line.decode("utf-8").split("\t")
                rsid = parts[1] + "_" + parts[0] #variant-gene
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