# eigensoft skript for all available NGS sequences
source("packages_load.R") #packages
source("01_sequences_for_alignment.R")

setwd("")
path_fasta = "./fastafiles_all/"
path_out = "./for_alignment/"

##Subtype B
sub_b_covar = preselect(subtype = "pooled", ART_naive = F, analysis = "HLA")
##read alignment
msa = Biostrings::readDNAMultipleAlignment("alignment_NT.fa",format = "fasta")
#Recode sequences
jk <- data.frame(id = names(unmasked(msa)), seq = unmasked(msa))
##Rename column names
colnames(jk) <-  c("id","seq")

##Dataset for sequencing data
maxfa = jk %>% separate(seq, into = paste("V", 1:(max(nchar(jk$seq))+1), sep = ""), sep = "")
fwrite(na="NA",maxfa,"maxfa.csv")
maxfa = maxfa %>%
  dplyr::select(-V1)
colnames(maxfa)[2:ncol(maxfa)] = paste0("V",1:(ncol(maxfa)-1))

#remove hxb2
maxfa = maxfa %>% dplyr::filter(!grepl(id, pattern = "HXB2|K03455.1")) %>% tibble()

##Merge sequencing data with patient data
gwas = merge(maxfa, sub_b_covar, by.x ="id", by.y = "base_uuid", all.x=FALSE, all.y = FALSE) %>% tibble()
colnames(gwas)[2:(ncol(maxfa))] = paste0("pos",seq(1,(ncol(maxfa)-1)))

##PLINK2
n_snps = length(maxfa)
n_samp = nrow(gwas)
vcf_allels = NULL
vcf_allels = as.data.frame(matrix(ncol = 1, nrow = n_snps-1))
for(i in 2:n_snps){
  vcf_allels$V1[i-1] = i-1
  vcf_allels$A[i-1] = sum(gwas[,i] == "A")
  vcf_allels$C[i-1] = sum(gwas[,i] == "C")
  vcf_allels$G[i-1] = sum(gwas[,i] == "G")
  vcf_allels$T[i-1] = sum(gwas[,i] == "T")
  vcf_allels$Z[i-1] = sum(gwas[,i] == "-")
}

##gwas as bases as characters
gwas_char = gwas
gwas_char[,2:n_snps][gwas_char[,2:n_snps] == "-"] = "Z"
fwrite(na="NA",gwas_char, "gwas_char_file.csv", row.names = FALSE)
##vcf file
gwas_vcf = NULL
gwas_vcf = as.data.frame(matrix(ncol = (n_samp+9), nrow = (n_snps-1)))
# maf should be 5 of all sequences
maf = 5
for(z in 1:(n_snps-1)){
  gwas_vcf[z,1] = 1
  gwas_vcf[z,2] = z
  gwas_vcf[z,3] = "."
  gwas_vcf[z,4] = colnames(vcf_allels)[which(vcf_allels[z,2:5] == max(vcf_allels[z,2:5]))[1]+1]
  gwas_vcf[z,5] = paste0(unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][1]),
                         ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][2]),
                         ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][3]),
                         ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][4]),
                         ",",unlist(colnames(vcf_allels)[which(vcf_allels[z,2:6] > maf)+1][5]))
  gwas_vcf[z,5] = gsub(pattern = paste0("NA,|NA|,NA"),"",gwas_vcf[z,5])
  gwas_vcf[z,5] = gsub(pattern = paste0("\\",gwas_vcf[z,4],","),"",gwas_vcf[z,5])
  gwas_vcf[z,5] = gsub(pattern = paste0(",\\",gwas_vcf[z,4]),"",gwas_vcf[z,5])
  gwas_vcf[z,5] = gsub(pattern = paste0("\\",gwas_vcf[z,4]),"",gwas_vcf[z,5])
  gwas_vcf[z,5][gwas_vcf[z,5] == ""] = "."
  gwas_vcf[z,6] = "."
  gwas_vcf[z,7] = "."
  gwas_vcf[z,8] = "."
  gwas_vcf[z,9] = "GT"

}
colnames(gwas_vcf)[1:9] = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

for(z in 1:n_samp){
  colnames(gwas_vcf)[z+9] = as.character(gwas_char$id[z])

  chrs = strsplit(as.character(paste0(gwas_vcf[1:(n_snps-1),4],",",
                                      gwas_vcf[1:(n_snps-1),5],",",
                                      gwas_char[z,2:n_snps])), split = ",")
  gwas_vcf[,z+9]  = paste0(unlist(lapply(chrs,anyDuplicated, fromLast = TRUE))-1)

  progress(z,progress.bar = FALSE, init = (z == 1), max.value = n_samp)
}

gwas_vcf[,][gwas_vcf[,] == "-1"] = "." # different from 02_binary_region_aa.R

###vcf binarisation
binarising = function(gwas_vcf,altc){
  othercs = c(1,2,3,4)[c(1,2,3,4) != altc]
  gwas_vcf_alt2 = gwas_vcf
  nseq = n_samp + 9
  gwas_vcf_alt2$ID = NA
  gwas_vcf_alt2$ID[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)] =
    paste0(gwas_vcf_alt2$POS[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)],
           gwas_vcf_alt2$REF[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc)],
           altc,sapply(strsplit(gwas_vcf_alt2[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc),5],","),"[[", altc))
  gwas_vcf_alt2 = gwas_vcf_alt2[which(!is.na(gwas_vcf_alt2$ID)),]
  gwas_vcf_alt2[,10:nseq][gwas_vcf_alt2[,10:nseq] == othercs[1]] ="."
  gwas_vcf_alt2[,10:nseq][gwas_vcf_alt2[,10:nseq] == othercs[2]] ="."
  gwas_vcf_alt2[,10:nseq][gwas_vcf_alt2[,10:nseq] == othercs[3]] ="."
  gwas_vcf_alt2[,10:nseq][gwas_vcf_alt2[,10:nseq] == altc] = 1
  gwas_vcf_alt2$ALT = sapply(strsplit(gwas_vcf_alt2[which(sapply(strsplit(gwas_vcf_alt2[1:(n_snps-1),5],","),length) >= altc),5],","),"[[", altc)
  return(gwas_vcf_alt2)
}
gwas_vcf_alt = rbind(binarising(gwas_vcf = gwas_vcf, altc = 1),
                     binarising(gwas_vcf = gwas_vcf, altc = 2),
                     binarising(gwas_vcf = gwas_vcf, altc = 3) ,
                     binarising(gwas_vcf = gwas_vcf, altc = 4))

## Exclude gaps, necessary for eigensoft
gwas_vcf_alt <- gwas_vcf_alt %>% dplyr::filter(ALT!="Z")

## needed for eigensoft
gwas_vcf_alt$POS <- c(1:nrow(gwas_vcf_alt))

fwrite(na="NA",gwas_vcf_alt,"plink2vcf.vcf", sep ="\t", quote = FALSE, row.names = FALSE)

header = paste0("##fileformat=VCFv1
##fileDate=",format(Sys.Date(),"%y%m%d"),"
##source=R",R.Version()$major,".",R.Version()$minor,"
##contig=<ID=1,length=15462>
##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
write.table(header,"vcfheader.txt",sep ="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##merge vcf with fileheader
system("cat vcfheader.txt  plink2vcf.vcf > vcfhead.vcf")

gwas_pheno = as.data.frame(sub_b_covar %>% ungroup() %>% dplyr::select(base_uuid, ART_naive)) 
colnames(gwas_pheno) = c("#IID","ART_naive")
gwas_pheno$ART_naive <- as.numeric(gwas_pheno$ART_naive)
fwrite(na="NA",gwas_pheno,"pheno_ART_naive.txt",sep ="\t", quote = FALSE, row.names = FALSE)

##GWAS
system("/Applications/plink2 --vcf vcfhead.vcf --chr-set -1 --make-pgen --sort-vars --out gwasp2")
system("/Applications/plink2 --pfile gwasp2 --pheno pheno_ART_naive.txt --nonfounders --make-pgen --out gwasp2")

system("/Applications/plink2 --pfile gwasp2 --mind 0.8 --nonfounders --make-pgen --out gwasp3") # 80% coverage
system("/Applications/plink2 --pfile gwasp3 --geno  0.05  --nonfounders --make-pgen --out gwasp4") # 90% non-missing

system("/Applications/plink2 --pfile gwasp4 --maf  0.05  --nonfounders --make-pgen --out gwasp5") # 10% minor allele frequency

## eigensoft
system("/Applications/plink2 --pfile gwasp5 --nonfounders --export ped --out eigensoft_gwas")
paste0("cd ", getwd())

#### run on terminal 
#### sed -i '' "s/-9/1/g" eigensoft_gwas.ped
#### conda activate eigensoft
#### convertf -p par.ped.eigenstrat.txt

pop_remove_eigensoft = fread("eigensoft_gwas.ind", sep= "", header = FALSE)
pop_remove_eigensoft$V3 = 1
fwrite(na="NA",pop_remove_eigensoft,"eigensoft_gwas.ind", sep= " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

#### smartpca -p parfile.txt

