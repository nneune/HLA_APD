# data preparation - the same for all three analyses ####
rm(list = ls())
source("packages_load.R") #packages
setwd("")
### loading the data set into R environment
...
# combine the data bases
...
# group the REGIONS
combined_db <- merge(combined_db, var_region, by="REGION") %>% dplyr::select(-REGION)
# all SHCS
all_SHCS <-  combined_db
# exclude those w/o VL measurements
combined_db <- combined_db %>% drop_na(RNA)
## ART-naive
# take first LABDATE as minlabdate
combined_db <- combined_db %>% group_by(ID) %>% dplyr::mutate(minlabdate = (min(LABDATE)))
# replacing HIV_POSDATE NA with HIV_POSDOCDATE
combined_db <- combined_db %>% 
  dplyr::mutate(HIV_POSDATE = ifelse(HIV_POSDATE >= HIV_POSDOCDATE | is.na(HIV_POSDATE), as.character(HIV_POSDOCDATE), as.character(HIV_POSDATE ))) %>% 
  dplyr::mutate(HIV_POSDATE=as.IDate(HIV_POSDATE, format =  "%Y-%m-%d")) 

# CD4 count over 300 & AIDS #
combined_db <- combined_db %>% dplyr::group_by(ID) %>% dplyr::mutate(AIDS_no = (minlabdate < FIRST_C_DATE -30 | is.na(FIRST_C_DATE)) & (CD4>300 | is.na(CD4))) # strict rules, no AIDS

# list of ART naive patients (min 30 days)
ART_ids <-combined_db %>% 
  filter(ART_START_DATE +7 >= LABDATE | is.na(ART_START_DATE)) %>%
  filter(minlabdate < FIRST_C_DATE -30 | is.na(FIRST_C_DATE)) %>% # no AIDS symptoms
  dplyr::summarise(n=n())

combined_db <- combined_db %>% dplyr::mutate(ART_naive = ID %in% ART_ids$ID) %>%
  relocate(minlabdate, .before = ART_START_DATE) %>%
  dplyr::select(-SEROC_DATE, -HIV_POSDOCDATE)

# NGS ####
NGS <- tibble(fread("")) 
NGS <- NGS %>% dplyr::select(ID=SHCS.ID,ZPHI.ID, sample.date,organism, base_uuid, PROJECT, sample.source) %>% dplyr::mutate(ZPHI.ID= ifelse(!grepl(ZPHI.ID,pattern="[[:digit:]]"),NA,ZPHI.ID))
# old NGS table compare if all data is present in new one
#transfer NA IDs from old to new, by ZPHI.ID
NGS <- full_join(
  NGS %>% distinct(),
  all_pot_uuids_220227 %>% dplyr::select(ID = SHCSID, ZPHI.ID = ZPHIID) %>% dplyr::mutate(ZPHI.ID = as.character(ZPHI.ID)) %>% drop_na(ZPHI.ID) %>% distinct(),
  by = "ZPHI.ID")  %>% 
  dplyr::mutate(ID = ifelse(is.na(ID.x), ID.y, ID.x)) %>% 
  dplyr::select(-ID.x,-ID.y) %>% 
  dplyr::arrange(ID)

NGS <- NGS %>% group_by(ID) %>% dplyr::arrange(ID) %>% dplyr::mutate(ZPHI.ID = as.character(mean(as.numeric(ZPHI.ID), na.rm=T))) %>% 
  dplyr::mutate(ZPHI.ID = ifelse(ZPHI.ID == "NaN",NA, ZPHI.ID))  

# transform sample date as date
NGS <- NGS %>% dplyr::mutate(sample.date = as.IDate(sample.date, format= "%d.%m.%Y")) 
# add subtypes (COMET and REGA) and APD scores (Carlise et al.) as references
NGS <- left_join(NGS, hiv_subtype %>% dplyr::select(subtype, base_uuid), by="base_uuid") %>% drop_na(ID) %>% dplyr::mutate(subtype = fct_infreq(subtype))
NGS <- left_join(NGS, APD_scores %>% dplyr::mutate(APD = ifelse(is.na(pol), gag,pol)) %>% dplyr::select(base_uuid, APD), by="base_uuid")

# clean up
NGS <- NGS %>%  #                                                          
  tidyr::drop_na(ID,sample.date) %>% # ID and sample date necessary       
  dplyr::filter(organism != "HCV") %>% # no Hepatitis C Virus Sequences   
  dplyr::distinct() %>% # no duplicated rows                        
  dplyr::filter(ID %in% all_SHCS$ID | ID %in% hlaBin$ID) %>% # real SHCS ID 
  dplyr::tibble()

# define ART_START_DATE and decide which sequences are ART naive
ART_dates <- distinct(all_SHCS %>% dplyr::select(ID, ART_START_DATE))

NGS <- merge(NGS,ART_dates,all=T) %>% distinct() %>% tibble()
NGS <- NGS %>% dplyr::mutate(ART_naive = sample.date <= ART_START_DATE+7 | is.na(ART_START_DATE)) %>% drop_na(base_uuid)

## filter out proviral seqs
NGS <- NGS %>% 
  dplyr::mutate_all(funs(replace(., .=='', NA))) %>%
  dplyr::mutate(sample.source = as.factor(sample.source),
                PROJECT = as.factor(PROJECT)) %>%
  dplyr::mutate(
    sample.source = fct_relevel(sample.source, "PBMC",after = Inf),
    sample.source = fct_relevel(sample.source, "cell",after = Inf),
    PROJECT = fct_relevel(PROJECT, "RetroSeq",after = Inf),
  )

NGS <- NGS %>%
  dplyr::mutate(proviral = ifelse(!is.na(sample.source) & (sample.source=="cell"|sample.source=="PBMC"),T,F)) %>% 
  dplyr::mutate(proviral = ifelse(grepl(PROJECT,pattern="ESS|LLV|Proviral|RetroSeq|Plasmid_sequencing|SingleGenomeProject_NFL"), T,proviral)) %>% 
  relocate(proviral, .before = PROJECT) 


# longitudinal data
NGS_ids_longi <- NGS %>% 
  filter(ART_naive) %>%
  group_by(ID) %>%  dplyr::summarise(n=n()) %>% filter(n>1)

NGS <- NGS %>% group_by(ID) %>% dplyr::mutate(in_longi = ID %in% NGS_ids_longi$ID)

NGS <- NGS %>% group_by(ID) %>%  dplyr::arrange(sample.date, .by_group = T) %>%
  dplyr::mutate(time_diff = difftime(sample.date, lag(sample.date), units = "days"))

# update my in_longi list
NGS_ids_longi <- NGS %>% filter(ART_naive) %>%group_by(ID) %>% dplyr::filter(mean(time_diff, na.rm=T) > 0 & in_longi ==T) %>%  dplyr::summarise(n=n())

NGS <- NGS %>% dplyr::mutate(in_longi = ID %in% NGS_ids_longi$ID)

# add viral PCs
NGS <- merge(NGS, eigensoft_gwas, by="base_uuid", all=T) %>% tibble()

# my original db all ids should include a yes/no statement
combined_db <- combined_db %>% group_by(ID) %>% dplyr::mutate(in_NGS = ID %in% NGS$ID)

combined_db <- combined_db %>% group_by(ID) %>%
  dplyr::mutate(in_NGS_longi = ID %in% NGS_ids_longi$ID) %>%
  dplyr::mutate(logRNA = log10(RNA+1)) # RNA log

# shipping over to ART-naive
chronic_ART_naive <- combined_db %>% filter(ART_naive==T) %>% dplyr::select(-ART_naive)

# actual SPVL
## mean and median
# time filter applied
SPVL_30 <- chronic_ART_naive %>% group_by(ID) %>%
  filter(LABDATE >  HIV_POSDATE + 30  |  LABDATE >  minlabdate +30 |  LABDATE >  REGDATE +30)%>%
  filter(LABDATE <= ART_START_DATE+7 | is.na(ART_START_DATE)) %>%
  filter(LABDATE < FIRST_C_DATE -30 | is.na(FIRST_C_DATE)) %>% # no AIDS symptoms
  dplyr::rename(REGION = VAR_DESC) %>%
  
  dplyr::mutate(mean_logRNA = mean(logRNA),
                median_logRNA = median(logRNA),
                median_CD4 = median(CD4, na.rm = T),
                sd_RNA= sd(logRNA),
                
                # redefining as factors
                REGION = as.factor(REGION),
                SEX = as.factor(SEX),
                ETHNICITY = as.factor(ETHNICITY),
                CENTER = as.factor(CENTER)) %>%
  
  dplyr::relocate(REGION, .before = ETHNICITY) %>%
  
  dplyr::mutate(time_since_infection = difftime(LABDATE,  HIV_POSDATE, units = "days")) %>%
  dplyr::mutate(time_since_infection =
                  ifelse((is.na(time_since_infection) & REGDATE + 30 <= minlabdate),
                         difftime(LABDATE,  REGDATE, units = "days"),
                         difftime(LABDATE,  HIV_POSDATE, units = "days"))) %>%
  
  dplyr::mutate(time_since_infection2 =
                  ifelse(is.na(time_since_infection),
                         difftime(LABDATE,  minlabdate, units = "days"),
                         time_since_infection)) %>%
  dplyr::select(-time_since_infection) %>%
  dplyr::rename(time_since_infection = time_since_infection2)

# strict rules
SPVL_30 <- SPVL_30 %>% 
  dplyr::mutate(strict = AIDS_no==T & ( sd_RNA<1.5 | is.na(sd_RNA)))

SPVL_30$REGION <- relevel(SPVL_30$REGION,ref = "Western Europe")


# sequences ####
maxf <- (fread("./ngs/csv_format_230508.csv", # generated in 00_formatting_threshhold_to_fasta.R
               stringsAsFactors=FALSE,colClasses = c("character")))[,-1] %>% tibble()
# all NAs are -
maxf[is.na(maxf)] <-  "-"

# collapse all columns
strings <- maxf %>% group_by(V1) %>%
  unite("fasta", V2: paste0("V",ncol(maxf)), sep = "")

# make fasta lines
strings <- strings %>% group_by(V1) %>%
  dplyr::mutate(comb = paste0(">", V1, "\n", fasta)) %>%
  dplyr::select(V1,comb)

colnames(strings) <- c("id", "fasta")

# how long are the sequences?
strings <- strings %>%  dplyr::select(base_uuid = id, fasta)
# 10'498 bases long
# 9'719 is the length of HXB2

# coverage
strings <- strings %>% group_by(base_uuid) %>%
  dplyr::mutate(coverage = sum(str_count(fasta, c("A", "T", "G", "C"))))

# filtering out ####
# filter earliest sample.date
# filter coverage if sample.date within 30 days
# filter > 0.1 coverage

# add NGS sample date and ID
NGS_s <- merge(strings, NGS, by = "base_uuid", all=T)

# filter 0.1
NGS_s <- NGS_s %>% group_by(ID) %>% 
  filter(coverage>=9719*0.1) # 9719 = length of hxb2

#blasting and codon alignment successful?
for (region in c("gag","pol","env","nef","vif","vpr","vpu", "rev", "tat")) {
  f_files <- list.files(paste0("./ngs/",region,"_codon_align"),pattern = "_AA.fa")
  assign(paste0("uuids_",region),unique(gsub(f_files,pattern=paste0("_",region,"_AA.fa"), replacement="")))
}

NGS_s <- NGS_s %>% filter(base_uuid %in% uuids_pol| 
                            base_uuid %in% uuids_nef| 
                            base_uuid %in% uuids_vif|
                            base_uuid %in% uuids_env|
                            base_uuid %in% uuids_gag|
                            base_uuid %in% uuids_rev|
                            base_uuid %in% uuids_tat|
                            base_uuid %in% uuids_vpr|
                            base_uuid %in% uuids_vpu)

NGS_s <- NGS_s %>% 
  dplyr::mutate(sample.date = as.Date(sample.date, format= "%Y-%m-%d")) %>%
  dplyr::mutate(ART_START_DATE = as.Date(ART_START_DATE, format= "%Y-%m-%d")) %>% 
  distinct()


# data dplyr::selection for alignment and hlaSNPbin ####
# choose ART naive yes or no, subtype all (no) or B (yes)
ART = c(T,F)
subtypes=c("B", "pooled") 

for(sub in subtypes){
  for (ART_naive in ART){
    cat(sub,ART_naive,"\n")
    
    # filter earliest
    NGS_strings <- NGS_s %>% 
      {if(sub == "nonB") filter(.,(subtype!="B")) else .} %>%
      {if(sub != "nonB" & sub != "pooled") filter(.,(subtype==sub)) else .} %>%
      {if(ART_naive == T) filter(.,(ART_naive==T)) else .} %>%
      group_by(ID) %>%
      relocate(ID, .before = base_uuid) %>% dplyr::select(-organism,-time_diff, -in_longi) %>% # longitudinal filter is taken out here - change for access to ids or go back to NGS data set
      # filter(min(sample.date) +30 >= sample.date) %>%   
      dplyr::arrange(desc(ART_naive), proviral,sample.date, desc(coverage), .by_group = T) %>% # highest coverage at earliest date, best if ART naive and not proviral DNA
      dplyr::slice(1) %>%
      ungroup() %>% tibble() %>%
      dplyr::mutate(subtype= fct_relevel(subtype,"B"))
    
    # NGS SPVL ####
    SPVL_strings <- merge(NGS_strings %>%
                            dplyr::select(-fasta,-ART_START_DATE),
                          SPVL_30 %>%
                            dplyr::select( -CD4, -RNA,-CD4DATE,-AIDS_no, -in_NGS), by="ID", all=T) %>% tibble()
    
    SPVL_strings <- distinct(SPVL_strings %>%
                               drop_na(base_uuid) %>% drop_na(mean_logRNA))
    
    SPVL_strings <- SPVL_strings  %>% group_by(ID) %>%
      dplyr:: dplyr::arrange(desc(coverage),desc(strict), desc(subtype), .by_group = T) %>%
      dplyr:: dplyr::arrange(LABDATE, .by_group = T)
    
    # generating phenotype for GWAS - can be used for all analysis containing NGS and VL
    # attention! SPVL and mean log ... are changed because of the 180 day rule - the LABDATE is only included if it is near the sample.date
    assign(paste(sep = "_","phenotype",sub,ifelse(ART_naive == T, "naive", "all")),
           SPVL_strings %>% 
             filter(ART_naive==T) %>%
             dplyr::filter(LABDATE >  HIV_POSDATE + 30  |  LABDATE >  minlabdate +30 |  LABDATE >  REGDATE +30) %>%
             dplyr::filter(LABDATE <= ART_START_DATE+7 | is.na(ART_START_DATE)) %>%
             dplyr::filter(LABDATE < FIRST_C_DATE -30 | is.na(FIRST_C_DATE)) %>% # no AIDS symptoms
             # ungroup() %>%
             dplyr::filter(LABDATE >= sample.date - 180 & LABDATE <= sample.date + 180) %>% # ensures that same virus strain is taken for SPVL
             # group_by(ID) %>%
             dplyr::mutate(mean_logRNA = mean(logRNA),
                           n_VL = n())  %>%
             dplyr::mutate(median_logRNA = median(logRNA)) %>%
             dplyr::mutate(sd_RNA= sd(logRNA)) %>%
             dplyr::mutate(subtype = as.factor(subtype),
                           SEX = as.factor(SEX),
                           ETHNICITY = as.factor(ETHNICITY),
                           CENTER = as.factor(CENTER),
                           REGION = as.factor(REGION)) %>%
             dplyr:: dplyr::arrange(desc(strict), desc(subtype), .by_group = T) %>%
             dplyr::slice(1) %>% # take the first
             ungroup()%>%
             dplyr::mutate(APDz = APD/ sd(APD, na.rm = T)))
    
    fwrite(x = get(paste(sep = "_","phenotype",sub,ifelse(ART_naive == T, "naive", "all"))), file = paste0("./phenotypes/SPVL_phenotype_",sub,ifelse(ART_naive == T, "_naive", "_all"),".csv"), na = "NA")
    
    # NGS HLA ####
    HLA_strings <- NGS_strings %>%
      dplyr::filter(ID %in% hlaBin$ID)
    
    HLA_strings <- HLA_strings  %>% group_by(ID) %>%
      dplyr:: dplyr::arrange(desc(coverage), desc(subtype), .by_group = T) %>%
      dplyr:: dplyr::arrange(sample.date, .by_group = T) %>%
      filter(sample.date == min(sample.date)) # keep only first sample date
    
    assign(paste0("hla_phenotype_",sub,ifelse(ART_naive == T, "_naive", "_all")),HLA_strings %>% group_by(ID) %>%
             dplyr:: dplyr::arrange(desc(coverage), desc(subtype), .by_group = T) %>%
             dplyr::slice(1) %>% dplyr::select(-coverage, -sample.date,-ART_START_DATE, -fasta)) 
    
    fwrite(x = get(paste0("hla_phenotype_",sub,ifelse(ART_naive == T, "_naive", "_all"))),file = paste0("./phenotypes/HLA_phenotype_",sub,ifelse(ART_naive == T, "_naive", "_all"),".csv"), na = "NA")
    
    
    # alignment ####
    path_fasta = "./ngs/"
    path_out = "./ngs/for_alignment"

    ##Subtype B
    sub_b_covar = preselect(subtype = sub, ART_naive = ART_naive,analysis = "HLA")

    ##AA
    for (region in c("gag","pol","env","nef","vif","vpr","vpu", "rev", "tat")){
    seq_foralign_selection(region = paste0(region),path_fasta = path_fasta,path_out = path_out,covar = sub_b_covar,type = "AA",subtype = sub,coverthres = 0.2,ART_naive = ART_naive,align_program ="mafft",longitudinal=  F)
      }
    
  }
}


## making gwas_aa_alt files ####
# sourcing file that converts the fasta files to a binary outcome with a csv containing all variants and all ids
generation_hlaSNPbin_file <- function(reg,keep_files=F,phenotype_binary, wd,sub,alignment_date = paste0(Sys.Date()), ART_naive = c(F,T)){
  if (!file.exists(paste0(wd))) {dir.create(paste0(wd))}
  setwd(paste0(wd)) 
  # running the two scripts to get binary file and hlaSNPbin for (un)dotted data
  alignment_date <<- alignment_date
  folder <<- paste0("./",reg,"/")
  cat(reg, " - ")
  f <- tempfile() # create a temp file
  template <- readLines("02_binary_region_aa.R") # read your template file
  writeLines(template, f) # write the template code into a temp connection
  gsub_file(f, "region", reg, fixed=TRUE) # rewrite the temp file
  gsub_file(f,']] ="."',paste0(']] = 0'), fixed=TRUE)
  gsub_file(f, '_aa_alt_dot', '_aa_alt', fixed=TRUE) # remove all dots, ref = all aa
  source(f) # run the temp file
  
  f <- tempfile() # create a temp file
  template <- readLines("03_hlaSNPbin_AA.R") # read your template file
  writeLines(template, f) # write the template code into a temp connection
  gsub_file(f, "region", reg, fixed=TRUE) # rewrite the temp file
  gsub_file(f, '_dot', '', fixed=TRUE)
  source(f) # run the temp file
  
}

subtypeB=F
keep_files = T # =FALSE: all files are deleted from environment
ART_naive = F

for(sub in subtypes){
  phenotype_binary <- get(paste0("hla_phenotype_",sub,ifelse(ART_naive == T, "_naive", "_all"))) %>% distinct()
  region_names= c("gag", "vif", "rev","vpu","vpr", "tat","nef", "pol","env" )
  for (region in region_names) {
    generation_hlaSNPbin_file(reg = region,phenotype_binary = phenotype_binary,alignment_date = Sys.Date(),ART_naive = F,sub = sub, wd = paste0("./ngs/hlaSNPbin/",sub,"/"))
  }
}

# importing all data in a fast way ####
## not for longitudinal data

# ART naive no is only for power analyses, but the hlaSNPbin should be ART naive anyway
ART_naive = T
subtypes="B"
for(sub in subtypes){
  cat("\n",sub, "\n")
  if(sub=="B"){F_test=T}else{F_test=F}
  setwd(dir = paste0("./ngs/hlaSNPbin/",sub,"/"))
  hxb2_ref <- tibble(fread("./hxb2_reference_nuclAA_gene.csv"))
  regions <- c("gag","pol","env","nef","vif","vpr", "rev", "tat","vpu")
  
  for (z in 1:length(regions)) {
    region = as.character(regions[z])
    cat("|",region)
    ## basics
    assign(paste0("hlaSNPbin_", region), tibble(fread(header = T,paste0("hlaSNPbin_", region,".csv"))))
    if(sum(grepl("pol",colnames(get(paste0("hlaSNPbin_", region)))))==0){
      assign(paste0("hlaSNPbin_", region), merge(get(paste0("hlaSNPbin_", region)), APD_scores, by = "base_uuid", all = T) %>% group_by(ID) %>% drop_na(ID))}
    
    # hlaSNPbin APD and PCs have to be standardized
    hlaSNPbin <- get(paste0("hlaSNPbin_",region)) %>% rowwise()%>%
      dplyr::mutate(APD = ifelse(region == "gag", gag, pol)) %>%
      ungroup() %>%
      dplyr::mutate(APDz = APD/sd(APD, na.rm = T)) %>%
      dplyr::select(-pol,-gag,-env) %>%
      dplyr::mutate(PC1z = PC1/sd(PC1, na.rm = T),
                    PC2z = PC2/sd(PC2, na.rm = T),
                    PC3z = PC3/sd(PC3, na.rm = T),
                    PC4z = PC4/sd(PC4, na.rm = T),
                    PC5z = PC5/sd(PC5, na.rm = T),
                    PC6z = PC6/sd(PC6, na.rm = T),
                    PC7z = PC7/sd(PC7, na.rm = T),
                    PC8z = PC8/sd(PC8, na.rm = T),
                    PC9z= PC9/sd(PC9, na.rm = T),
                    PC10z = PC10/sd(PC10, na.rm = T))
    
    # make a translation key
    translation <-  tibble(fread(paste0("./",region,"/SNP_",region,"_table_ref.csv")))[2:4]
    gwasalt_hxb2_ref <- tibble(fread(paste0(region,"/seqs_to_keep_",region,".csv")))
    translation <- merge(translation, gwasalt_hxb2_ref %>% dplyr::select(-Var2), by.x="POS",by.y="gwas_pos", all=T)
    translation$insertion = NA
    translation$multins = ".0"
    translation$index = NA
    if(is.na(translation[1,5])){
      translation[1,5] <- 0
      translation[1,6] <- "0_1_ins"
      translation[1,7] <- ".1"
    }
    for (i in 1:nrow(translation)) { # this loop creates information if and where inserts occur
      if (is.na(translation$hxb2_index[i])) {
        translation$index[i] = translation$index[i-1]
        translation$insertion[i] = paste0(translation$index[i],"_",translation$index[i]+1,"_ins")
        translation$multins[i] = case_when(
          translation$index[i] == translation$index[i-4] & translation$POS[i] != translation$POS[i-1]~ paste0(".",str_sub(translation$multins[i-1], start=2) %>% as.numeric() + 1),
          translation$index[i] == translation$index[i-4] & translation$POS[i] == translation$POS[i-1]~ translation$multins[i-1],
          translation$index[i] == translation$index[i-3] & translation$POS[i] != translation$POS[i-1]~ paste0(".",str_sub(translation$multins[i-1], start=2) %>% as.numeric() + 1),
          translation$index[i] == translation$index[i-3] & translation$POS[i] == translation$POS[i-1]~ translation$multins[i-1],
          translation$index[i] == translation$index[i-2] & translation$POS[i] != translation$POS[i-1]~ paste0(".",str_sub(translation$multins[i-1], start=2) %>% as.numeric() + 1),
          translation$index[i] == translation$index[i-2] & translation$POS[i] == translation$POS[i-1]~ ".1",
          translation$index[i] == translation$index[i-1] ~ paste0(".1"))
      }else{translation$index[i] =translation$hxb2_index[i]}
    }
    translation <- translation %>% dplyr::mutate(multins = ifelse(multins == ".0", NA, multins))
    last(colnames(translation)) <- "region_no"
    colnames(hxb2_ref)[2] <- "reg"
    translation <- merge(translation,hxb2_ref %>% filter(reg == region), by="region_no")
    colnames(translation) <- c("protein_no_hbx2", "POS","ID","id","hxb2_ref","hxb2_index","insertion","multins","hxb2_pos_nucl","region","gene","gene_no")
    translation <- translation %>% drop_na(ID) %>% tibble()
    translation <- translation %>%
      dplyr::mutate(var_name_long = paste0(toupper(region), protein_no_hbx2,
                                           gsub(id,pattern="1|2|3|4|5|6|7|8|9|0", replacement = "")),
                    REF_ALT = str_remove_all(id,pattern="1|2|3|4|5|6|7|8|9|0")) %>%
      rowwise()%>%
      dplyr::mutate(REF = unlist(str_split(REF_ALT, ""))[1],
                    ALT = unlist(str_split(REF_ALT, ""))[2]) %>% dplyr::select(-REF_ALT) %>%
      dplyr::mutate(var_name_short = paste0(REF,protein_no_hbx2,ALT),
                    name_expl = ifelse(is.na(insertion), 
                                       paste0("x",protein_no_hbx2, ALT), 
                                       paste0("x",insertion,multins,".",ALT))) %>% ungroup()
    
    # replace the names in the hlaSNPbin with the real names
    hlaSNPbin_named <- hlaSNPbin
    
    #change names of hlaSNPbin
    col_indices <- match(colnames(hlaSNPbin_named), translation$ID)
    colnames(hlaSNPbin_named)[!is.na(col_indices)] <- translation$name_expl[na.omit(col_indices)]
    fwrite(hlaSNPbin_named, paste0("hlaSNPbin_named_", region,".csv"), row.names = F, na = "NA")
    
    #change names of Power-Fisher Table
    ## power analysis import
    if(F_test==T){
      assign(paste0("Fisher_", region), tibble(fread(paste0("1_Fisher_pvalues_corrected_",region,"_0.8_0.2.csv"))))
      Fisher_named <- get(paste0("Fisher_", region))
      col_indices <- match(colnames(Fisher_named), translation$ID)
      colnames(Fisher_named)[!is.na(col_indices)] <- translation$name_expl[na.omit(col_indices)]
      fwrite(Fisher_named, paste0("Fisher_named_", region,".csv"), row.names = F, na = "NA")
      rm("Fisher_named")
    }
    
    rm(list=c("translation", "hlaSNPbin","hlaSNPbin_named"))
  }
  # remove unspecified data frames, can lead to problems in function after
}

## How many combinations are left after Power analysis and Fisher tests + correction?
cat(sum(!is.na(Fisher_gag[,-1]),!is.na(Fisher_pol[,-1]),!is.na(Fisher_vif[,-1]),!is.na(Fisher_vpr[,-1]),!is.na(Fisher_tat[,-1]),!is.na(Fisher_rev[,-1]),!is.na(Fisher_vpu[,-1]),!is.na(Fisher_env[,-1]),!is.na(Fisher_nef[,-1])))

# sign interactions ####
# loop over all Fisher, then do glms for each row in the table and select those with sign. interaction -> that should be the new list of sign_hits!!
for (sub in subtypes){
  df = data.frame(NA,NA,NA,NA)
  colnames(df) <- (c("hla_allele", "var"  , "gene", "POS"))
  ART_naive=T
  count=0
  excluded_var <- NULL
  # Fisher from all
  # hlaSNPbin from ART-Naive
  setwd("")
  for (region in c("vif", "vpr", "tat","rev","vpu","nef", "gag", "env", "pol")) {
    Fisher <- fread(paste0("./ngs/hlaSNPbin/B/Fisher_named_",region,".csv")) %>% tibble()
    if (nrow(Fisher)<1) {next}
    sign_hits <- Fisher %>% pivot_longer(cols = !hla_allele, names_to = "var", values_to = "FDR") %>% drop_na(FDR) 
    hlaSNPbin <- fread(paste0("./ngs/hlaSNPbin/",sub,"/hlaSNPbin_named_",region,".csv")) %>% tibble()
    hlaSNPbin <- left_join(hlaSNPbin %>% distinct(), NGS_s %>%ungroup()%>% dplyr::select(base_uuid,subtype,ART_naive, starts_with("vPC"))) %>% tibble()
    hlaSNPbin$subtype <- relevel(as.factor(hlaSNPbin$subtype), ref=ifelse(sub=="B", sub, as.factor(hlaSNPbin$subtype) %>% table() %>% sort() %>% last() %>% names()))
    translation <- fread(paste0("./ngs/hlaSNPbin/",sub,"/translation_named_", region,".csv")) %>% tibble()
    
    # filter out high APD and remove sequences which are not ART-naive
    hlaSNPbin <- hlaSNPbin %>% filter(APD<0.05) %>% tibble() # filter APD<0.05 because of non-linear effect of APD on VL
    hlaSNPbin <- hlaSNPbin %>% filter(ART_naive==T) %>% tibble() # filter if they are ART naive
    hlaSNPbin <- hlaSNPbin %>%  ungroup() %>%
      dplyr::mutate(APDz = APD/sd(APD, na.rm = T)) %>%
      dplyr::mutate(PC1z = PC1/sd(PC1, na.rm = T),
                    PC2z = PC2/sd(PC2, na.rm = T),
                    PC3z = PC3/sd(PC3, na.rm = T),
                    PC4z = PC4/sd(PC4, na.rm = T),
                    PC5z = PC5/sd(PC5, na.rm = T),
                    PC6z = PC6/sd(PC6, na.rm = T),
                    PC7z = PC7/sd(PC7, na.rm = T),
                    PC8z = PC8/sd(PC8, na.rm = T),
                    PC9z= PC9/sd(PC9, na.rm = T),
                    PC10z = PC10/sd(PC10, na.rm = T)) %>%      
      dplyr::mutate(vPC1z = vPC1/sd(vPC1, na.rm = T), # this only for nonB subtypes
                    vPC2z = vPC2/sd(vPC2, na.rm = T),
                    vPC3z = vPC3/sd(vPC3, na.rm = T),
                    vPC4z = vPC4/sd(vPC4, na.rm = T),
                    vPC5z = vPC5/sd(vPC5, na.rm = T),
                    vPC6z = vPC6/sd(vPC6, na.rm = T),
                    vPC7z = vPC7/sd(vPC7, na.rm = T),
                    vPC8z = vPC8/sd(vPC8, na.rm = T),
                    vPC9z= vPC9/sd(vPC9, na.rm = T),
                    vPC10z = vPC10/sd(vPC10, na.rm = T))
    
    count=nrow(sign_hits %>% filter(!var %in% colnames(hlaSNPbin)))+count
    excluded_var <- rbind(excluded_var,sign_hits %>% filter(!var %in% colnames(hlaSNPbin)) %>% dplyr::mutate(gene = region))
    
    sign_hits <- sign_hits %>% filter(var %in% colnames(hlaSNPbin))
    sign_hits <- sign_hits %>% dplyr::mutate_if(is.character, as.factor)
    
    sign_hits <- sign_hits %>%  dplyr::arrange(FDR) %>% rowwise() %>%
      dplyr::mutate(coef = list(glm(data= hlaSNPbin, 
                                    formula= paste0(var,"~",hla_allele, "*APDz+PC1z+PC2z+PC3z+PC4z+PC5z+PC6z+PC7z+PC8z+PC9z+PC10z", 
                                                    ifelse(sub!="B","+vPC1z+vPC2z+vPC3z+vPC4z+vPC5z+vPC6z+vPC7z+vPC8z+vPC9z+vPC10z","")),
                                    family = "binomial") 
                                %>% tidy())) %>% dplyr::mutate(pval = (coef[["p.value"]] %>% last())) %>%
      dplyr::filter(pval <=0.05) %>% dplyr::select(-coef) %>% ungroup() %>% dplyr::mutate(gene = region) # significant interaction p-value (<0.05) is filtered
    
    sign_hits <- merge(sign_hits,translation %>% dplyr::select(var=name_expl,POS=hxb2_pos_nucl)) %>%  dplyr::arrange(POS)
    
    assign(paste0("hlaSNPvl_", region),full_join(get(paste0("phenotype_",sub,"_",ifelse(ART_naive == T, "naive", "all"))) %>% ungroup() %>% dplyr::select(-ID,-APD,-APDz),hlaSNPbin %>% ungroup() %>% dplyr::select(-ID), by="base_uuid") %>%
             drop_na(mean_logRNA,PC1))
    
    assign(paste0("hlaSNPbin_", region), hlaSNPbin)
    df <- full_join(df, sign_hits)
  }
  assign(paste0("count_",sub),count) # how many not found, should be == 0
  sign_hits <- df%>%  dplyr::arrange(POS) %>% drop_na() %>% dplyr::rename(hla = hla_allele)
  # all significant interactions with APD are selected
  fwrite(sign_hits, paste0("./phenotypes/sign_int_",sub,"_hits_",ifelse(ART_naive == T, "naive", "all"),".csv"), row.names = F, na = "NA")
  assign(paste0("sign_hits_",sub), fread(paste0("./phenotypes/sign_int_",sub,"_hits_",ifelse(ART_naive == T, "naive", "all"),".csv")))
}

# longitudinal ####
# IDs in HLA and NGS longi
longi_NGS_HLA <- merge(NGS_s%>% dplyr::select(-fasta,-time_diff) %>% distinct() ,hlaBin %>%  distinct(), by=c("ID"), all=T) %>% 
  filter(ART_naive ==T) %>% drop_na(sample.date, PC1, PC2) %>% distinct() %>% tibble()

sub = "pooled"
longi_NGS_HLA <- longi_NGS_HLA %>% drop_na(coverage) %>% filter(coverage>9719*0.1) %>% # exclude those without fasta & low coverage
  group_by(ID) %>%  dplyr::arrange(desc(coverage),.by_group = T)

longi_NGS_HLA$subtype <- relevel(as.factor(longi_NGS_HLA$subtype), ref="B")
longi_NGS_HLA <- longi_NGS_HLA %>% 
  {if(sub != "pooled") filter(.,ifelse(sub=="B",subtype == "B",subtype != "B")) else . } %>%
  filter(ART_naive ==T) #%>% filter(proviral ==F)

longi_NGS_HLA <- longi_NGS_HLA %>%  dplyr::arrange(desc(coverage))%>% group_by(ID) %>% distinct(sample.date, .keep_all = T)
longi_NGS_HLA <- longi_NGS_HLA %>% dplyr::filter(n()>1) %>% dplyr::select(-in_longi)

longi_NGS_HLA <- longi_NGS_HLA %>% group_by(ID)%>% dplyr::arrange(sample.date,.by_group = T) %>% 
  dplyr::mutate(time_lag = difftime(sample.date, lag(sample.date), units = "days")) %>% 
  dplyr::mutate(time = difftime(sample.date, min(sample.date), units = "days")) %>% 
  relocate(time, time_lag, .after=sample.date) %>% 
  dplyr::mutate(time = ifelse(is.na(time), 0, time))

fwrite(paste0("./ngs/longitudinal/longi_NGS_HLA_",sub,"_naive.csv"),x = longi_NGS_HLA, row.names = F, na = "NA")


longi_NGS_HLA <- longi_NGS_HLA %>% filter(any(time_lag>=30)) # sequences must be 30 days apart

# meta dataframe ####
meta_combined_db <- tibble(full_join(all_SHCS, NGS %>% dplyr::select(-ART_START_DATE), by="ID",relationship = "many-to-many")) %>% 
  drop_na(base_uuid)%>% 
  group_by(ID) %>% 
  dplyr::mutate(minlabdate = (min(LABDATE))) %>% 
  dplyr::mutate(HIV_POSDATE = ifelse(HIV_POSDATE >= HIV_POSDOCDATE | is.na(HIV_POSDATE), as.character(HIV_POSDOCDATE), as.character(HIV_POSDATE ))) %>% 
  dplyr::mutate(HIV_POSDATE=as.IDate(HIV_POSDATE, format =  "%Y-%m-%d"))  %>%
  dplyr::select(-contains("PC"), -RNA,-LABDATE,-CD4DATE,-FIRST_C_DATE,-SEROC_DATE,-CD4,-HIV_POSDOCDATE,-HIV_NEGDATE) %>% ungroup() %>% distinct()

meta_combined_db <- tibble(full_join(phenotype_pooled_all %>% dplyr::select(-contains("PC"),-ZPHI.ID,-proviral,-PROJECT,-sample.source,-subtype,-APD,-ART_naive,-contains("DATE")) %>% dplyr::mutate_if(is.Date, as.IDate),
                                     distinct(meta_combined_db) %>% dplyr::mutate_if(is.Date, as.IDate) %>% 
                                       dplyr::mutate(
                                         REGION = relevel(as.factor(VAR_DESC),ref="Western Europe"),
                                         ETHNICITY = relevel(as.factor(ETHNICITY),ref="1"),
                                         SEX = relevel(as.factor(SEX),ref="1"),
                                         RISKGROUP = relevel(as.factor(RISKGROUP),ref="MSM"),
                                         CENTER = relevel(as.factor(CENTER),ref="10")))) %>% distinct()%>%  dplyr::arrange(ID,base_uuid) %>% dplyr::select(-VAR_DESC)

meta_combined_db <-meta_combined_db %>%
  dplyr::mutate(VL = ifelse(is.na(mean_logRNA),F,T), # VL measurement
                HLA = ifelse(ID %in% hlaBin$ID, T, F), # HLA data
                NGS = ifelse(base_uuid %in% hla_phenotype_pooled_all$base_uuid, T, F), # True for all (phenotype dplyr::selects for it)
                longi = ifelse(ID %in% longi_NGS_HLA$ID, T, F)) %>%
  group_by(ID) %>%
  dplyr::mutate(date_diagnosis = as.Date(min(HIV_POSDATE,minlabdate,REGDATE,ART_START_DATE, sample.date, na.rm = T))) %>%
  relocate(HIV_POSDATE,REGDATE,minlabdate,ART_START_DATE, sample.date, .before = date_diagnosis)

meta_combined_db <- meta_combined_db %>% tibble() %>% ungroup() %>% dplyr::mutate(subtype=relevel(as.factor(subtype),ref="B")) %>%
  dplyr::mutate(ETHNICITY=set_labels(ETHNICITY, labels = c("other","White", "Black", "Hispano-american","Asian","other"))) %>% distinct() 

meta_combined_db <- 
  meta_combined_db %>%
  dplyr::mutate(`Baseline Age` = year(REGDATE) -BORN) %>%
  dplyr::mutate(REGION =fct_infreq(REGION)) %>% 
  dplyr::mutate(subtype = ifelse(subtype == "B", "B", "non-B")) %>%
  dplyr::mutate(subtype = relevel(as.factor(subtype), "B")) %>%
  dplyr::mutate(`standardized APD` = APD/sd(APD,na.rm = T)) %>%
  dplyr::mutate(`Risk group` = fct_infreq(RISKGROUP)) %>%
  dplyr::rename(
    `mean log10 RNA` = mean_logRNA,
    `median log10 RNA`=median_logRNA,
    `median CD4` = median_CD4,
    `ART naive` = ART_naive) %>%
  distinct() %>% ungroup() %>% dplyr::select(-RISKGROUP) 

## PERINAT have wrong "diagnosis" date
meta_combined_db <- meta_combined_db %>% 
  group_by(ID) %>%
  dplyr::mutate(date_diagnosis = as.Date(ifelse(is.na(HIV_POSDATE) & `Risk group`=="PERINAT", paste0(BORN, "-01-01"),date_diagnosis)))

# pairs of interest ####
pair_hits <- (get(paste0("sign_hits_B")) %>%  dplyr::mutate(NAME = paste0(hla,var)) %>% 
                filter(
                  NAME=="B_5701x33R"|
                    NAME=="A_0301x432R"|
                    NAME=="C_0501x223V"|
                    NAME=="B_4001x57E"|
                    NAME=="B_5101x37V"|
                    NAME=="A_1101x85V"|
                    NAME=="A_3201x704L"|
                    NAME=="B_4402x726D"|
                    NAME=="A_2402x135F"))[,1:4] %>% rbind(., data.frame(
                      hla="B_5701",
                      var="x242N",
                      gene="gag",
                      POS=1513
                    )) %>%  dplyr::arrange(gene)

pair_hits$gene <- factor(pair_hits$gene, levels = c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))

pair_hits <- pair_hits %>%  dplyr::arrange(gene) %>%
  dplyr::mutate(col= c("#332288","#0066CC", "#12A4EA","#117733","#EFB412","#C00000","#CC0066", "#731941"))


