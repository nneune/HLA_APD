#### run on cluster
folder <- getwd()
username = Sys.getenv()['USER'] # necessary to call packages
library(dplyr, lib.loc=paste("/data/",username,"/rpackages",sep=""))
library(readr, lib.loc=paste("/data/",username,"/rpackages",sep=""))
library(vcd, lib.loc=paste("/data/",username,"/rpackages",sep=""))
library(qdapTools, lib.loc=paste("/data/",username,"/rpackages",sep=""))
library(Hmisc, lib.loc=paste("/data/",username,"/rpackages",sep=""))
library(tibble, lib.loc=paste("/data/",username,"/rpackages",sep=""))
library(foreach, lib.loc=paste("/data/",username,"/rpackages",sep=""))
library(doParallel, lib.loc=paste("/data/",username,"/rpackages",sep=""))
library(data.table, lib.loc=paste("/data/",username,"/rpackages",sep=""))
#power calculation
get.power<-function(x, y){
  l=table(x,y)
  as.numeric(Hmisc::bpower(p1= (l[3]/(l[3] + l[1])),p2= (l[4]/(l[4] + l[2])),
                           n1=(l[3] + l[1]),
                           n2=(l[4] + l[2]),
                           odds.ratio = OR,
                           alpha=pval))}
# Fisher's exact test
get.fisher<-function(x, y){
  l=table(x,y)
  test <- stats::fisher.test(x=as.matrix(l))
  as.numeric(test$p.value)
}
cat("\nPackages and Functions loaded\n")

# coefficients 
pval = 0.05; cat("alpha = ", pval)
OR = 3; cat("\nOR = ", OR)
power_th = 0.8; cat("\npower threshold = ", power_th)
FDR_th = 0.2; cat("\nFDR threshold = ", FDR_th,"\n")
# loop
region_names= c("vif", "vpr", "tat","rev","nef", "gag", "pol","vpu","env")
registerDoParallel(cores = 9) # a core per gene # check if supported
foreach(region = region_names)%dopar%{
  region = as.character(region)
  hlaSNPbin = fread(paste0(folder,"/hlaSNPbin_",region,".csv")) %>% tibble()
  cat("\nLoading worked -", region,"\n")
  hlaSNPbin[3:ncol(hlaSNPbin)] <- mutate_if(hlaSNPbin[3:ncol(hlaSNPbin)], is.character, as.numeric)
  
  # power ####
  SNPs <- hlaSNPbin[grep("V",names(hlaSNPbin))]
  HLAs <- hlaSNPbin[grep("_",names(hlaSNPbin))][,-1]
  cat("Before:: Number of SNPs: ",ncol(SNPs), "  Number of HLA alleles: ", ncol(HLAs), "\n", " Number of theoretical combinations: ", ncol(HLAs)* ncol(SNPs),"\n")
  SNPs <- SNPs %>%
    select_if(colSums(., na.rm=T)>0.005*nrow(hlaSNPbin) & colSums(., na.rm=T)<nrow(hlaSNPbin)) %>% # minimal variant frequency of 0.5%
    select_if(colMeans(is.na(.)) < 0.2) # less than 20% missing
  HLAs <- HLAs %>% select_if(colSums(., na.rm=T)>1 & colSums(., na.rm=T)<nrow(hlaSNPbin))
  d = as.matrix(SNPs)
  r = as.matrix(HLAs)
  d.cols <- split(d, col(d))
  r.cols <- split(r, col(r))
  power_matrix = outer(r.cols,d.cols, Vectorize(get.power))
  power_matrix <- as.data.frame(power_matrix)
  rownames(power_matrix)<- names(HLAs)
  colnames(power_matrix)<- names(SNPs)
  cat(paste0("Power matrix finished - ", region,"\n"))
  fwrite(na="NA",power_matrix,file = paste0(folder,"/power_analysis_",region,"_05_",OR,".csv"),row.names = T)
  
  # Fisher ####
  power_matrix = fread(paste0(folder,"/power_analysis_",region,"_05_",OR,".csv"))
  Fisher_pvalues <- NULL
  names(power_matrix)[1] <- "hla_allele"
  power_matrix <- as.data.frame(power_matrix)
  power_matrix <- power_matrix %>% group_by(hla_allele) %>% dplyr::replace(power_matrix <= power_th, NA) # remove all value below power threshold (0.8)
  cat("Before:  Number of HLA alleles: ", nrow(power_matrix), "  Number of variants: ", ncol(power_matrix), "  Region: ",region,"\n")
  # HLA excluded which NA
  power_matrix <- power_matrix[which((rowSums(is.na(power_matrix[,-1])) == ncol(power_matrix[,-1])) ==F),]
  # variants excluded which NA
  power_matrix <- power_matrix[,-(which(colSums(is.na(power_matrix[,-1])) == nrow(power_matrix))+1)]
  SNPs <- hlaSNPbin[grep("V",names(hlaSNPbin))]
  HLAs <- hlaSNPbin[grep("_",names(hlaSNPbin))][,-1]
  SNPs <- SNPs %>% select(colnames(power_matrix[grepl("V",names(power_matrix) )]))
  HLAs <- HLAs %>% select(power_matrix$hla_allele)
  cat("Number of SNPs: ",ncol(SNPs), "  Number of HLA alleles: ", ncol(HLAs), "\n", " Number of combinations: ", ncol(HLAs)* ncol(SNPs),"\n")
  d = as.matrix(SNPs)
  r = as.matrix(HLAs)
  d.cols <- split(d, col(d))
  r.cols <- split(r, col(r))
  Fisher_pvalues = outer(r.cols,d.cols, Vectorize(get.fisher))
  Fisher_pvalues <- as.data.frame(Fisher_pvalues)
  rownames(Fisher_pvalues)<- names(HLAs)
  colnames(Fisher_pvalues)<- names(SNPs)
  cat("Fisher finished - ",region,"\n")
  fwrite(na="NA",Fisher_pvalues,file = paste0(folder,"/Fisher_pvalues_",region,"_",power_th,"_",OR,".csv"),row.names = T)
  
  ## power analysis import
  powertable <- power_matrix
  cat(sum(!is.na(powertable[,-1])), "pairs survived in", region, ": power\n") # how many variants
  powertable[,-1] <- replace(powertable[,-1], is.na(powertable[,-1]), 0) # setting all to 0 and 1s for matrix multiplication
  powertable[,-1] <- replace(powertable[,-1], powertable[,-1]!=0, 1)
  
  # Fisher multiple testing ####
  Fisher_pvalues <- as.data.frame(Fisher_pvalues) %>% rownames_to_column("hla_allele") %>% tibble()
  Fisher_matrix <- as.matrix(powertable %>% column_to_rownames("hla_allele"))*as.matrix(Fisher_pvalues%>% column_to_rownames("hla_allele")) # matrix multiplication
  Fisher_matrix <- Fisher_matrix %>% replace(Fisher_matrix ==0, NA)
  Fisher_pvalues<- as.data.frame(Fisher_matrix) %>% rownames_to_column("hla_allele") %>% group_by(hla_allele)
  cat(sum(!is.na(Fisher_pvalues %>% ungroup() %>% dplyr::select(-hla_allele))), "pairs survived in", region, ": Fischer\n") # how many variants
  ## P value correction "Benjamini Hochberg" aka FDR
  # Filter the columns with names containing "V"
  filtered_cols <- grepl("V", names(Fisher_pvalues))
  Fisher_cols <- Fisher_pvalues[, filtered_cols]
  # Apply p.adjust() function to the filtered columns
  adjusted_cols <- matrix(p.adjust(as.matrix(Fisher_cols), method = "BH", n = sum(!is.na(Fisher_pvalues[, -1]))),ncol=ncol(Fisher_cols))
  # Replace the filtered columns with the adjusted values
  Fisher_pvalues[, filtered_cols] <- adjusted_cols
  Fisher_pvalues<- as.data.frame(Fisher_pvalues) %>%  column_to_rownames("hla_allele")
  Fisher_pvalues<- Fisher_pvalues%>% dplyr::replace(Fisher_pvalues> FDR_th, NA) %>% rownames_to_column("hla_allele") # FDR p-value significance threshold to 0.2
  if(nrow(Fisher_pvalues)>1){
    cat(sum(!is.na(Fisher_pvalues[,-1])), "pairs survived in", region, ": FDR\n") # how many variants
    Fisher_pvalues<- Fisher_pvalues[which((rowSums(is.na(Fisher_pvalues[,-1])) == ncol(Fisher_pvalues[,-1])) ==F),]
    Fisher_pvalues<- Fisher_pvalues[,-(which(colSums(is.na(Fisher_pvalues[,-1])) == nrow(Fisher_pvalues))+1)]
    fwrite(na="NA",Fisher_pvalues,file = paste0(folder,"/Fisher_pvalues_corrected_",region,"_",power_th,"_",FDR_th,".csv"),row.names = F)
  }else{print(sum(!is.na(Fisher_pvalues[,-1]))," pairs survived FDR - ", region)}
  cat(region, "END\n",format(Sys.time(), format = "%H:%M:%S"),"\n________________________________________________________________________________________________________________________________________________________________________")
}
