setwd("")

##libraries
library(devtools)
library(tidyverse)
library(foreach)
library(doParallel)

devtools::install_github("mhahsler/rBLAST", force = TRUE)
library(rBLAST)

##functions file
source("01_apd_functions.R")

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "./ncbi-blast-2.15.0+/bin/", sep= .Platform$path.sep))
Sys.which("blastn")

##sequence selection
##required: frequency file, coverage file, consensus sequence
##specify path to these files
freq_path =  "./freq/"
cov_path =  "./cov/"
consensus_path ="./fastafiles_all/"

##specify path to MACSE (codon alignment)
macse_path = "./macse/macse_v2.05.jar"

##define list of uuids
uuids <- str_remove_all(list.files(consensus_path), pattern="\\.fa$")

uuids = as.data.frame(uuids)
colnames(uuids) = "uuid"

uuids <- drop_na(uuids)
uuids <- remove_rownames(uuids) %>% distinct()


## reference HXB2; nucleotide positions taken for analysis
# `5'LTR` = c(1, 634)
# gag = c(790, 2292)
# pol = c(2085, 5096)
# vif = c(5041, 5619)
# vpr = c(5559, 5850)
# tat = rbind(c(5831, 6045), c(8379, 8469))
# rev = rbind(c(5970, 6045), c(8379, 8653))
# vpu = c(6062,6310)
# env= c(6225,8795)
# nef= c(8797,9417)
# `3'LTR` = c(9086, 9719)
region_names <- c("gag","pol","env","nef","vif","vpr","vpu", 
                  "rev_1", "rev_2", "tat_1","tat_2") # rev and tat are split

##extract regions via blast (reference hxb2)
registerDoParallel(cores = 4)
foreach(region = region_names)%dopar%{ 
  region_extract(uuids$uuid,region)}

##codon alignment of blast vs reference (may take a lot of time with several sequences)
registerDoParallel(cores = 4)
foreach(region = region_names)%dopar%{ 
  condon_align(uuids$uuid,region)}

# successful blast uuids needed
uuids_pol = fread(paste0("./blastmeta/blastmeta_pol.csv")) %>% filter(!is.na(SubjectID)) %>% select(uuid = QueryID)
uuids_gag = fread(paste0("./blastmeta/blastmeta_gag.csv")) %>% filter(!is.na(SubjectID)) %>% select(uuid = QueryID)
uuids_env = fread(paste0("./blastmeta/blastmeta_env.csv")) %>% filter(!is.na(SubjectID)) %>% select(uuid = QueryID)

##calculated apd score; for genes pol, env and gag (see Carlisle et al, 2019)
apd_scores = uuids %>% 
  rowwise() %>%
  dplyr::mutate(pol = ifelse((uuid %in% uuids_pol$uuid), 
                             apd_calc(uuid = uuid,region = "pol", depth_threshold = 20), # or 100, depending on threshold preference
                             NA)) %>%
  dplyr::mutate(env = ifelse((uuid %in% uuids_env$uuid), 
                             apd_calc(uuid = uuid,region = "env", depth_threshold = 20), 
                             NA)) %>%
  dplyr::mutate(gag = ifelse((uuid %in% uuids_gag$uuid),
                             apd_calc(uuid = uuid,region = "gag", depth_threshold = 20),
                             NA))

write.csv(apd_scores,"./apd_scores_all.csv")

