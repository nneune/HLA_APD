##power AA bin ## called within 11_data_load.R ## region is replaced by respective HIV gene
# loading RNA AA alignment
source("packages_load.R") #packages
cat("region")
if (ART_naive==T) {
  cat("ART naive\n")
  msa_region = readAAMultipleAlignment(paste0("./alignment/mafft/nci_AA_",sub,"_mafft_region_ART_aligned_",alignment_date,".fa"))
}
if (ART_naive==F) {
  cat("all sequences\n")
  msa_region = readAAMultipleAlignment(paste0("./alignment/mafft/nci_AA_",sub,"_mafft_region_all_aligned_",alignment_date,".fa"))
}
cat("Sequence loaded\n")
#Recode sequences
jk <- data.frame(id = names(unmasked(msa_region)), seq = unmasked(msa_region))
##Rename column names
colnames(jk) <-  c("id","seq")

if (file.exists(folder)) {cat("\n")} else {dir.create(folder)}

##Dataset for sequencing data
maxfa = jk %>% separate(seq, into = paste("V", 1:(max(nchar(jk$seq))+1), sep = ""), sep = "")
fwrite(na="NA",maxfa,paste0(folder,"maxfa_region.csv"))
maxfa = maxfa %>%
  dplyr::select(-V1)
colnames(maxfa)[2:ncol(maxfa)] = paste0("V",1:(ncol(maxfa)-1))

#remove non HXB2 positions
seqs_to_keep_region = as.matrix(maxfa %>% dplyr::filter(grepl(id, pattern = "HXB2")))
seqs_to_keep_region  = reshape2::melt(seqs_to_keep_region)

seqs_to_keep_region  = seqs_to_keep_region  %>%
  dplyr::select(-Var1) %>%
  dplyr::filter(!grepl(value, pattern = "HXB2")) %>%
  dplyr::filter(value != "-") %>%
  dplyr::mutate(hxb2_index_region = as.numeric(rownames(.)))

seqs_to_keep_region$gwas_pos = as.numeric(gsub("V", "", seqs_to_keep_region$Var2))
fwrite(na="NA",seqs_to_keep_region, paste0(folder,"seqs_to_keep_region.csv"))
cat("translation: seqs_to_keep_region written to ", getwd(),"\n")
#remove hxb2
maxfa = maxfa %>% dplyr::filter(!grepl(id, pattern = "HXB2"))

##Merge sequencing data with patient data
# phenotype_binary <- preselect("B")
gwas = merge(maxfa, phenotype_binary, by.x ="id", by.y = "base_uuid", all.x=FALSE, all.y = FALSE) # phenotype_binary is defined by Master script
colnames(gwas)[2:(ncol(maxfa))] = paste0("pos",seq(1,(ncol(maxfa)-1)))

##dplyr::filter subtype B
cat("subtypedplyr::filtered",ifelse(subtypeB==T,"", "not"),"for 'B':",dim(gwas%>% group_by(ID)%>%filter(!subtype == "B"))[1],"of",dim(gwas)[1],"removed\n")
gwas <- gwas %>% group_by(ID) %>% {if (subtypeB==T)dplyr::filter(.,subtype == "B") else . }

cat("Kept from alignment by phenotype:",dim(gwas)[1],"of",nrow(phenotype_binary),"\n")

vcf_allels_region_aa = NULL
vcf_allels_region_aa = as.data.frame(matrix(ncol = 1, nrow = (ncol(maxfa)-1)))

gwas_region_aa = gwas
for(i in 2:(ncol(maxfa))){
  vcf_allels_region_aa$V1[i-1] = i - 1
  vcf_allels_region_aa$F[i-1] = sum(gwas_region_aa[,i] == "F")
  vcf_allels_region_aa$L[i-1] = sum(gwas_region_aa[,i] == "L")
  vcf_allels_region_aa$I[i-1] = sum(gwas_region_aa[,i] == "I")
  vcf_allels_region_aa$M[i-1] = sum(gwas_region_aa[,i] == "M")
  vcf_allels_region_aa$V[i-1] = sum(gwas_region_aa[,i] == "V")
  vcf_allels_region_aa$S[i-1] = sum(gwas_region_aa[,i] == "S")
  vcf_allels_region_aa$P[i-1] = sum(gwas_region_aa[,i] == "P")
  vcf_allels_region_aa$T[i-1] = sum(gwas_region_aa[,i] == "T")
  vcf_allels_region_aa$A[i-1] = sum(gwas_region_aa[,i] == "A")
  vcf_allels_region_aa$Y[i-1] = sum(gwas_region_aa[,i] == "Y")
  vcf_allels_region_aa$H[i-1] = sum(gwas_region_aa[,i] == "H")
  vcf_allels_region_aa$Q[i-1] = sum(gwas_region_aa[,i] == "Q")
  vcf_allels_region_aa$N[i-1] = sum(gwas_region_aa[,i] == "N")
  vcf_allels_region_aa$K[i-1] = sum(gwas_region_aa[,i] == "K")
  vcf_allels_region_aa$D[i-1] = sum(gwas_region_aa[,i] == "D")
  vcf_allels_region_aa$E[i-1] = sum(gwas_region_aa[,i] == "E")
  vcf_allels_region_aa$C[i-1] = sum(gwas_region_aa[,i] == "C")
  vcf_allels_region_aa$W[i-1] = sum(gwas_region_aa[,i] == "W")
  vcf_allels_region_aa$R[i-1] = sum(gwas_region_aa[,i] == "R")
  vcf_allels_region_aa$G[i-1] = sum(gwas_region_aa[,i] == "G")
  vcf_allels_region_aa$Z[i-1] = sum(gwas_region_aa[,i] == "-")
  vcf_allels_region_aa$X[i-1] = sum(gwas_region_aa[,i] == "X")
}
vcf_allels_region_aa$var = apply(vcf_allels_region_aa[2:22],1, function(x) sum(x>1))

##gwas as bases as characters
gwas_region_char_aa = gwas_region_aa[,c(1,2:ncol(maxfa))]
gwas_region_char_aa[,2:ncol(maxfa)][gwas_region_char_aa[,2:ncol(maxfa)] == "-"] = NA

##vcf file
gwas_region_vcf_aa = NULL
gwas_region_vcf_aa = as.data.frame(matrix(ncol = nrow(gwas)+9, nrow = (ncol(maxfa)-1)))
maf = 5 # major allele frequency threshold
for(z in 1:(ncol(maxfa)-1)){
  gwas_region_vcf_aa[z,1] = 1
  gwas_region_vcf_aa[z,2] = z
  gwas_region_vcf_aa[z,3] = "."
  gwas_region_vcf_aa[z,4] = colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:21] == max(vcf_allels_region_aa[z,2:21]))[1]+1]
  gwas_region_vcf_aa[z,5] = paste0(unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][1]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf )+1][2]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][3]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][4]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][5]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][6]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][7]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][8]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][9]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][10]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][11]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][12]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][13]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][14]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][15]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][16]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][17]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][18]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][19]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][20]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][21]),
                                   ",",unlist(colnames(vcf_allels_region_aa)[which(vcf_allels_region_aa[z,2:22] > maf)+1][22]))
  gwas_region_vcf_aa[z,5] = gsub(pattern = paste0("NA,|NA|,NA"),"",gwas_region_vcf_aa[z,5])
  gwas_region_vcf_aa[z,5] = gsub(pattern = paste0(gwas_region_vcf_aa[z,4],","),"",gwas_region_vcf_aa[z,5])
  gwas_region_vcf_aa[z,5] = gsub(pattern = paste0(",",gwas_region_vcf_aa[z,4]),"",gwas_region_vcf_aa[z,5])
  gwas_region_vcf_aa[z,5] = gsub(pattern = paste0(gwas_region_vcf_aa[z,4]),"",gwas_region_vcf_aa[z,5])
  gwas_region_vcf_aa[z,5][gwas_region_vcf_aa[z,5] == ""] = NA
  gwas_region_vcf_aa[z,6] = "."
  gwas_region_vcf_aa[z,7] = "."
  gwas_region_vcf_aa[z,8] = "."
  gwas_region_vcf_aa[z,9] = "GT"
  
}
colnames(gwas_region_vcf_aa)[1:9] = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

library(svMisc)
for(z in 1:length(gwas_region_aa$id)){
  colnames(gwas_region_vcf_aa)[z+9] = as.character(gwas_region_char_aa[z,1])
  chrs = strsplit(as.character(paste0(gwas_region_vcf_aa[1:(ncol(maxfa)-1),4],",",
                                      gwas_region_vcf_aa[1:(ncol(maxfa)-1),5],",",
                                      gwas_region_char_aa[z,2:(ncol(maxfa))])), split = ",")
  gwas_region_vcf_aa[,z+9]  = paste0(unlist(lapply(chrs,anyDuplicated, fromLast = TRUE))-1)
  
  progress(z,progress.bar = T, init = (z == 1), max.value = length(gwas_region_aa$id))
}

gwas_region_vcf_aa[,][gwas_region_vcf_aa[,] == "-1"] = NA

###vcf binarisation
binarising = function(gwas_region_vcf_aa,altc){
  nseq = nrow(gwas) + 9
  othercs = c(1,2,3,4,5,6,7,8,9,10,11,12,
              13,14,15,16,17,18,19,20,21,22)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22) != altc]
  gwas_region_aa_alt2 = gwas_region_vcf_aa
  gwas_region_aa_alt2$ID = NA
  gwas_region_aa_alt2$ID[which(sapply(strsplit(gwas_region_aa_alt2[1:(ncol(maxfa)-1),5],","),length) >= altc)] =
    paste0(gwas_region_aa_alt2$POS[which(sapply(strsplit(gwas_region_aa_alt2[1:(ncol(maxfa)-1),5],","),length) >= altc)],
           gwas_region_aa_alt2$REF[which(sapply(strsplit(gwas_region_aa_alt2[1:(ncol(maxfa)-1),5],","),length) >= altc)], # REF takes the consensus not HXB2!!
           altc,sapply(strsplit(gwas_region_aa_alt2[which(sapply(strsplit(gwas_region_aa_alt2[1:(ncol(maxfa)-1),5],","),length) >= altc),5],","),"[[", altc))
  gwas_region_aa_alt2 = gwas_region_aa_alt2[which(!is.na(gwas_region_aa_alt2$ID)),]
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[1]] ="."    # "."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[2]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[3]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[4]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[5]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[6]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[7]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[8]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[9]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[10]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[11]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[12]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[13]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[14]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[15]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[16]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[17]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[18]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[19]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[20]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[21]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == othercs[22]] ="."
  gwas_region_aa_alt2[,10:nseq][gwas_region_aa_alt2[,10:nseq] == altc] = 1
  gwas_region_aa_alt2$ALT = sapply(strsplit(gwas_region_aa_alt2[which(
    sapply(strsplit(gwas_region_aa_alt2[1:(ncol(maxfa)-1),5],","),length) >= altc),5],","),"[[", altc)
  return(gwas_region_aa_alt2)
}
binary_region_aa_alt = NULL
binary_region_aa_alt = rbind(binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 1),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 2),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 3),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 4),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 5),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 6),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 7),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 8),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 9),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 10),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 11),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 12),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 13),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 14),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 15),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 16),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 17),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 18),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 19),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 20),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 21),
                             binarising(gwas_region_vcf_aa = gwas_region_vcf_aa, altc = 22))

binary_region_aa_alt <- binary_region_aa_alt %>% drop_na(ALT) %>% select(-QUAL, -INFO, -FILTER, -`#CHROM`, -FORMAT) %>% arrange(POS)
fwrite(na="NA",x = binary_region_aa_alt,file = paste0(folder,"binary_region_aa_alt_dot"), sep ="\t", quote = FALSE, row.names = FALSE)
cat("\nsaved to:", getwd(),folder)
rm("z", "i","jk")

if(keep_files == F){rm("maxfa", "maf", 
                       "gwas_region_vcf_aa","gwas_region_char_aa","gwas","vcf_allels_region_aa", "msa_region","gwas_region_aa")
}

cat("\n_______________________________________________________________________________________________________________________________________________________\n")
