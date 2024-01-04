rm(list = ls())

#Loading necessary libraries
source("packages_load.R") #packages
library(sqldf)
library(io)
# BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
# library(VCFArray)
library(fs)
library(qqman)
# library(GWASTools)
library(Rfast)
library(devtools)
library(magrittr)
library(nnet)

setwd("")

##Loading csv file from Sandras pipeline which has the probability for each base at each position
##Take the base with the highest probability
##List all files 
frequenciescsv_files= vector()
frequenciescsv_files = list.files("var_thresh_files")
##generate the file pathways
frequenciespaths = vector()
frequenciespaths = paste0("var_thresh_files/", frequenciescsv_files)

# folder paths
freq_path="./freq_files/"
cov_path="./cov_files/"

##Extract the uuids of each file as a identifier for each sequence later
uuids=frequenciescsv_files %>% str_sub(start=ifelse(nchar(.)>80,nchar("base_uuid= "),1),end= ifelse(nchar(.)>80,nchar(.)-nchar(".merged_ITERATION_COUNT_variant_threshold.csv"),nchar(.)))

##Dataframe with uuids and pathway
presel = as.data.frame(cbind(uuids,frequenciespaths))
presel <- presel[!uuids %in% (list.files("fastafiles_all/") %>% str_sub(end=nchar(.)-3)),] # those which are missing

why = 0 
##warnings not relevant
for(z in 1:length(presel$uuids)){
  if(presel$uuids[z] %in% str_remove(list.files(freq_path), "_freq.csv")){next} # if already present, skip
  marytesat = fread(file=presel$frequenciespaths[z])
  if(is.na(colnames(marytesat)[2]) || colnames(marytesat)[2] != "POS"){ 
    why = why + 1
    next
  }
  marytesat = marytesat %>%
    separate(ALT, c("first", "second","third","fourth"), sep = "/", remove = FALSE) %>%
    separate(AF, c("first_f", "second_f","third_f","fourth_f"), sep = "/", remove = FALSE) %>%
    mutate(first_f = as.numeric(first_f),second_f = as.numeric(second_f),third_f = as.numeric(third_f),fourth_f = as.numeric(fourth_f)) %>%
    mutate(ref_f = 100 - rowSums(dplyr::select(., first_f, second_f, third_f,fourth_f), na.rm = TRUE)) %>%
    mutate(A = case_when(REF == "A" ~ paste0(ref_f),
                         first == "A" ~ paste0(first_f), 
                         second == "A" ~ paste0(second_f),
                         third == "A" ~ paste0(third_f),
                         fourth == "A" ~ paste0(fourth_f))) %>%
    mutate(T = case_when(REF == "T" ~ paste0(ref_f),
                         first == "T" ~ paste0(first_f), 
                         second == "T" ~ paste0(second_f),
                         third == "T" ~ paste0(third_f),
                         fourth == "T" ~ paste0(fourth_f))) %>%
    mutate(G = case_when(REF == "G" ~ paste0(ref_f),
                         first == "G" ~ paste0(first_f), 
                         second == "G" ~ paste0(second_f),
                         third == "G" ~ paste0(third_f),
                         fourth == "G" ~ paste0(fourth_f))) %>%
    mutate(C = case_when(REF == "C" ~ paste0(ref_f),
                         first == "C" ~ paste0(first_f), 
                         second == "C" ~ paste0(second_f),
                         third == "C" ~ paste0(third_f),
                         fourth == "C" ~ paste0(fourth_f))) %>%
    mutate(GAP = NA) %>%
    dplyr::select(POS,A,C,G,T,GAP,COV)
  
  marytesat = marytesat %>%
    melt(., id.vars = c("POS","COV")) %>%
    group_by(POS) %>%
    arrange(match(variable, c("A", "C", "G","T","GAP")),.by_group = TRUE) 
  
  marytesat$variable = as.character(marytesat$variable)
  marytesat$variable[which(marytesat$variable == "GAP")] = "-"
  marytesat$value[which(is.na(marytesat$value))] = 0
  marytesat$value = as.numeric(marytesat$value) / 100
  
  freq <- marytesat %>%
    dplyr::select(POS,variable,value)
  
  colnames(freq) = c("position_1basedindexing", "variant", "frequency")
  fwrite(freq,paste0(freq_path,presel$uuids[z],"_freq.csv"), row.names = FALSE, na="NA")
  
  cov <- marytesat %>%
    mutate(reference = "K03455.1") %>%
    dplyr::select(reference,POS,COV) %>%
    distinct(reference,POS,COV)
  
  colnames(cov) = c("reference", "position", "depth")
  fwrite(cov,paste0(cov_path,presel$uuids[z],"_cov.csv"), row.names = FALSE, na="NA")
  
  progress(z,progress.bar = F, init = (z == 1), max.value = length(presel$uuids))
}

##List all files 
frequenciescsv_files= vector()
frequenciescsv_files = list.files(freq_path)
##generate the file pathways
frequenciespaths = vector()
frequenciespaths = paste0(freq_path, frequenciescsv_files)
##Get filenames of coverage files
frequenciescsv_files_av= vector()
frequenciescsv_files_av = list.files(cov_path)
##Generate datapaths of coverage files
frequenciespaths_av = vector()
frequenciespaths_av = paste0(cov_path, frequenciescsv_files_av)

##Extract the uuids of each file as a identifier for each sequence later
uuids = str_remove(list.files(freq_path), "_freq.csv")
uuids_av = str_remove(list.files(cov_path), "_cov.csv")
##Dataframe with uuids and pathway
presel = as.data.frame(cbind(uuids,uuids_av,frequenciespaths,frequenciespaths_av))


##Defining data frame for sequences (9719 basepairs + identifier)
maxf = as.data.frame(matrix(ncol = 10500, nrow = length(presel$uuids)))
##coding all variables as -> missing
##Write sequence identifier in sequence data frame
maxf[,1] = presel$uuids
##if index instead of base is returned: (- = 1 A = 2 C = 3 G = 4 T = 5)
##no implementation to return bases under threshold yet
##Main loop
for(z in 1:length(presel$uuids)){
  ##load frequency file
  cte = fread(file = presel$frequenciespaths[z])
  ##load sequence depth file
  avg = fread(file = presel$frequenciespaths_av[z])  
  ##change from long to wide format
  ctw_wide <- spread(cte, position_1basedindexing, frequency)
  ##Cache data frame
  hxb2pos = data.frame(pos=as.numeric(array(1:10499))) 
  ##calculates max frequency and returns the base, then merges with Base position !!Potential initial time lack, does not occur when the index instead of the Base is returned (index to base can be done later)
  hxb2pos = full_join(hxb2pos,
                  data.frame(pos=as.numeric(colnames(ctw_wide)[2:length(ctw_wide)]),
                             V2=ctw_wide[apply(ctw_wide[,2:length(ctw_wide)],2,nnet::which.is.max),1]),
                  by = "pos")
  ##merge cache with the average sequence depth
  hxb2pos = full_join(hxb2pos,avg %>% dplyr::rename(pos=position), by="pos")
  ##actual position variable after filtering for sequencing depth
  ##apply sequence depth filter
  hxb2pos <- hxb2pos %>% dplyr::mutate(acpos = ifelse(depth > 20, V2, NA))
  #write into data frame
  maxf[z,2:10500] = hxb2pos$acpos
  if(!is.na(maxf[z,"V41"])){
    if(maxf[z,"V41"] == TRUE){
      print(z)
      print(maxf$V1[z])
      break
    }
  }
  progress(z,progress.bar = FALSE, init = (z == 1), 
           max.value = length(presel$uuids))
}


##Save sequences as csv
fwrite(maxf,paste0("csv_format_",format(Sys.Date(),"%y%m%d"),".csv"),na="NA")


# sequences ####
##might be some error when reading in, some positions are "TRUE"
maxf = fread(paste0("csv_format_",format(Sys.Date(),"%y%m%d"),".csv"), stringsAsFactors=FALSE, 
                  colClasses = c("character"))

# all NAs are -
maxf[is.na(maxf)] = "-"

# collapse all columns
strings <- maxf %>% group_by(V1) %>%
  unite("fasta", V2: paste0("V",ncol(maxf)), sep = "")

# make fasta lines
strings <- strings %>% group_by(V1) %>%
  dplyr::mutate(comb = paste0(">", V1, "\n", fasta)) %>%
  dplyr::select(V1,comb)

# make fasta files
strings %>% filter(V1 %!in% (list.files("fastafiles_all") %>% str_sub(end=nchar(.)-3))) %>% group_by(V1) %>%
  group_walk(~write.table(.x,paste0("fastafiles_all/",.y$V1, ".fa"),col.names = F, quote = F, row.names = F))

system("cat ./fastafiles_all/*.fa > ./fastafiles_all/merged_fasta_all_.fa")
system("sed '/^[^>]/s/[N]/-/g' ./fastafiles_all/merged_fasta_all_.fa > ./fastafiles_all/merged_fasta_all.fa")
