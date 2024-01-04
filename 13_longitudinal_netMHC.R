## loading packages and data
setwd("")
source("packages_load.R")

# alignment ####
path_fasta = "./fasta/"
path_out = "./for_alignment/"

## data loaded from 11_data_load.R
phenotype_binary = longi_NGS_HLA %>% filter(ART_naive,subtype=="B") %>% group_by(ID) %>% filter(n()>1)
region_names <- c("gag","pol","env","nef","vif","vpr","vpu", "rev", "tat")

for (region in region_names){
  seq_foralign_selection(paste0(region),path_fasta,path_out,covar = phenotype_binary,type = "AA",subtype = "B",
                         coverthres = 0.2,ART_naive = T,align_program ="mafft", longitudinal = T)
}


# hlaSNPbin generation ####
registerDoParallel(cores = 9)
foreach(region = region_names)%dopar%{
  generation_hlaSNPbin_file(region,phenotype_binary,alignment_date = Sys.Date(),ART_naive = T, wd = "./hlaSNPbin/")
}

# data prep for survival analysis ####
ART_naive = T
setwd("")
hxb2_ref <- tibble(fread("hxb2_reference_nuclAA_gene.csv")) 
F_test=F
# different dataset - longitudinal 
for (z in 1:length(region_names)) {
  region = as.character(region_names[z])
  cat("|",region)
  ## basics
  hlaSNPbin <-  tibble(fread(header = T,paste0("hlaSNPbin_", region,".csv")))
  if(sum(grepl("pol",colnames(hlaSNPbin)))==0){
    hlaSNPbin <-  merge(hlaSNPbin, APD_scores, by = "base_uuid", all = T) %>%
      group_by(ID) %>% drop_na(ID)}
  
  hlaSNPbin <- hlaSNPbin %>% rowwise()%>%
    dplyr::mutate(APD = ifelse(region == "gag", gag, pol)) %>%
    ungroup() %>%
    dplyr::select(-pol,-gag,-env) 
  
  hlaSNPbin <- right_join(NGS_s %>%ungroup()%>% dplyr::select(base_uuid,subtype,sample.date,ART_naive, proviral),
                          hlaSNPbin %>% distinct(), 
                          by = join_by(base_uuid)) %>% tibble()
  hlaSNPbin$subtype <- fct_infreq(hlaSNPbin$subtype)  
  # make a translation key
  translation <-    tibble(fread(paste0("./",region,"/SNP_",region,"_table_ref.csv"))[,2:4])
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
  for (i in 1:nrow(translation)) {
    
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
  translation <- translation %>% mutate(multins = ifelse(multins == ".0", NA, multins))
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
                                     paste0("x",insertion,multins,".",ALT))
    ) %>% ungroup()
  
  # assign(paste0("translation_", region), translation)
  fwrite(translation, paste0("translation_named_", region,".csv"), row.names = F, na = "NA")
  
  # replace the names in the hlaSNPbin with the real names
  hlaSNPbin_named <- hlaSNPbin
  
  #change names of hlaSNPbin
  col_indices <- match(colnames(hlaSNPbin_named), translation$ID)
  colnames(hlaSNPbin_named)[!is.na(col_indices)] <- translation$name_expl[na.omit(col_indices)]
  
  fwrite(na = "NA",hlaSNPbin_named, paste0("hlaSNPbin_named_longi_", region,".csv"), row.names = F)
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
}


# Analysis part ####
sign_hits_all <- (fread("hlaSNPbin/model_APD_OR_B.csv") %>% tibble())[1:3] %>% distinct() 

# remove single IDs
for (region in region_names) {
  assign(paste0("hlaSNPbin_longi_", region), fread(paste0("./hlaSNPbin/hlaSNPbin_named_longi_", region,".csv")) %>% 
           filter(subtype=="B" & ART_naive) %>% select(-ART_naive) %>%
           group_by(ID) %>% 
           filter(n()>1) %>% 
           arrange(sample.date,.by_group = T) %>% 
           mutate(time = difftime(sample.date, min(sample.date), units = "days")) %>% 
           relocate(time,.after=sample.date) %>% 
           mutate(time = ifelse(is.na(time), 0, time))
  )
  
}


table_var_change <- sign_hits_all %>% rowwise() %>% 
  dplyr::mutate(var = list(get(paste0("hlaSNPbin_longi_", gene)) %>% dplyr::select(ID, sample.date, time,var=any_of(VAR),hla=any_of(HLA)) %>% 
                             group_by(ID) %>%
                             dplyr::arrange((sample.date), .by_group = T) %>% 
                             drop_na(var,hla) %>% # remove NA
                             dplyr::mutate(time = difftime(sample.date, min(sample.date), units = "days")) %>% 
                             dplyr::filter(any(var == 0 & time==0)) %>%
                             dplyr::filter(any(var == 1 & time!=0)) %>%
                             dplyr::filter(time!=0 & var!=0) %>% 
                             dplyr::filter(time==min(time)) %>% 
                             ungroup() %>%
                             mutate(time=as.numeric(time/365.25)) %>%
                             dplyr::summarise(time.median= median(time),
                                              time.sum= sum(time),
                                              time.max= max(time), n=n()))) %>% 
  
  dplyr::mutate(hla = list(get(paste0("hlaSNPbin_longi_", gene)) %>% dplyr::select(ID, sample.date, var=any_of(VAR),hla=any_of(HLA), time) %>%
                             group_by(ID) %>%
                             dplyr::arrange((sample.date), .by_group = T) %>% 
                             drop_na(var,hla) %>% # remove NA
                             dplyr::mutate(time = difftime(sample.date, min(sample.date), units = "days")) %>% 
                             dplyr::filter(any(var == 0 & time==0)) %>%
                             # dplyr::filter(any(var == 1 & time!=0)) %>%
                             dplyr::filter(n()>1) %>%
                             dplyr::filter(time!=0) %>% 
                             # dplyr::filter(time==min(time)) %>% 
                             dplyr::mutate(
                               time1 = ifelse(var == 1, time, max(time)),
                               time1 = min(time1)) %>%
                             dplyr::filter(time == time1) %>%
                             mutate(time=as.numeric(time/365.25)) %>%
                             ungroup() %>%
                             group_by(hla) %>%
                             dplyr::summarise(time.median= median(time),
                                              time.sum= sum(time),
                                              var.sum = sum(var),
                                              time.max= max(time),
                                              n=n(),
                                              inci.rate= var.sum/time.sum*100))) %>%
  unnest_wider(var, names_sep = ".") %>% 
  unnest_wider(hla, names_sep = ".") %>%
  unnest() %>%
  dplyr::mutate(do_longi= var.n>2)  # filter that at least 3 events (changes from 0 to 1) happen (independent of HLA allele)


fwrite(na = "NA",table_var_change,"./table_var_change.csv",row.names = F)


setwd("")
#### survival analysis function ####
# a function to get the KM & netMHCpan results with the variant and HLA input
longitudinal_analysis <- function(VAR,HLA,GENE,mosaic_time,color="black",pair=F,longi=F, plotting="off",DO_LONGI = F){
  gene_variant <- paste0(firstup(GENE), str_remove(VAR, pattern = "x"))
  variant_label <- str_remove(VAR, pattern = "x")
  hla_label <- insert_colon(str_replace_all(HLA, "_", "*"))
  
  # netMHCpan
  netMHCpan_rank <- NULL
  if(nchar(HLA)<7){ # only class I
    # dplyr::selecting netMHCpan 9-mer peptide
    seqs_to_keep_hxb2 <- fread(paste0("./hlaSNPbin/",GENE,"/seqs_to_keep_",GENE,".csv")) %>% rowwise() %>%
      mutate(Var2 = paste0("V",as.numeric(str_split(Var2,pattern="V")[[1]][2])+1)) %>% ungroup()
    
    maxfa <- fread(paste0("./hlaSNPbin/",GENE,"/maxfa_",GENE,".csv"))%>% tibble()
    pos <- fread(paste0("./hlaSNPbin/translation_named_", GENE,".csv")) %>% filter(name_expl == VAR) %$% protein_no_hbx2
    
    # maxfa seq position
    vars <- seqs_to_keep_hxb2[seqs_to_keep_hxb2 %$% get(paste0("hxb2_index_",GENE)) %in% c((pos-8):(pos+8)),]$Var2
    
    hxb2_cons <- str_c(collapse = "", unlist(maxfa[vars][1,]))
    
    consensus <- apply(maxfa[vars], 2, function(x) {
      freq_table <- table(x)
      max_freq <- max(freq_table)
      letters_with_max_freq <- names(freq_table)[freq_table == max_freq]
      return(letters_with_max_freq)
    })
    
    string_cons <- str_c(collapse = "", unlist(consensus))
    consensus[9] = str_remove_all(VAR, pattern = "x[[:digit:]]+")
    string_mut <- str_c(collapse = "", unlist(consensus))
    
    x = paste0(">Consensus\n", string_cons,"\n>Mutation\n", string_mut)
    folder = paste0("./netMHCpan_analyses/sign_hits")
    f = paste0("peptide_",gene_variant,"_",HLA,".txt")
    
    if (!file.exists(folder)) {dir.create(folder)}
    write_file(x,file = paste0(folder,"/",f))
    
    log_file <- system(paste0("cd ",folder," \n       ../netMHCpan -f ",f," -BA -xls -a HLA-",insert_colon(str_remove(HLA, pattern="_"))," -l 9 "),intern=T)
    writeLines(log_file, paste0(folder,"/NetMHCpan_",gene_variant,"_",HLA,"_out.txt"))
    
    suppressWarnings(log_file <- read_table(paste0(folder,"/NetMHCpan_",gene_variant,"_",HLA,"_out.txt"), 
                                            skip = 51, comment = "##",show_col_types = F))
    
    log_file$Identity <- c(rep("Consensus",nrow(log_file)/2), rep("Mutation",nrow(log_file)/2))
    log_file= log_file %>% drop_na(MHC) %>% filter(grepl(Pos,pattern="[[:digit:]]"))  
    log_file[12:16] = mutate_if(log_file[12:16], is.character, as.numeric)
    log_file= log_file %>% rowwise() %>% 
      mutate(
        strong = ifelse(`%Rank_EL` <= 0.500,1,0),
        weak = ifelse(`%Rank_EL` <= 2.000 && strong != 1,1,0)) %>% 
      dplyr::select(-BindLevel)
    
    fwrite(na="NA",log_file, paste0(folder,"/NetMHCpan_",gene_variant,"_",HLA,"_out.txt"))
    
    log_file$Pos <- c((nrow(log_file)/2):1,(nrow(log_file)/2):1)
    
    netMHCpan <- log_file %>% ungroup() %>%
      dplyr::select(Pos, MHC,ID=Identity, Peptide,EL_Rank= `%Rank_EL`, 
                    Aff = `Aff(nM)`, strong , weak) %>% group_by(Pos) %>%
      filter(mean(strong) != 0 | mean(weak) != 0) %>%
      mutate(Pos = paste("position", Pos))
    
    if(nrow(netMHCpan)!=0){
      netMHCpan <- netMHCpan %>% arrange(Pos, ID)
      
      netMHCpan[c("EL_Rank","Aff")] <- mutate_if(netMHCpan[c("EL_Rank","Aff")], is.character, as.numeric)
      
      
      a <- ggplot(data=netMHCpan, aes(y=Aff, x=Pos, fill=ID)) +
        geom_bar(position=position_dodge(width=0.8), stat="identity", color="black",width = 0.8)+
        scale_fill_manual("Variant", values= c(colorspace::lighten(color, 0.6),color))+
        theme_classic()+
        xlab(label= "Position in 9-mer") +
        ylab(label = "MHC Binding Affinity (nM)")+
        theme(
          axis.ticks.length = unit(2, "mm"),
          axis.ticks.y = element_line(color="black"),
          axis.text = element_text(size=14),
          title = element_text(size=14, face = "bold"),
          axis.title = element_text(size=15, face = "bold"))+
        ggtitle(paste0(gene_variant, ":", hla_label))
      
      
      b <- ggplot(data=netMHCpan, aes(y=EL_Rank, x=Pos, fill=ID)) +
        geom_bar(position=position_dodge(width=0.8), stat="identity", color="black",width = 0.8)+
        scale_fill_manual("Variant",values= c(colorspace::lighten(color, 0.6),color))+
        {if(pair==T) guides(fill = guide_legend(override.aes = list(fill = c("grey80","grey40"))))}+
        theme_classic()+
        scale_y_continuous(breaks=c(0,0.5,2,5,7.5,10, 12.5,15))+
        xlab(label= "Position in 9-mer") +
        ylab(label = "Mass Spectrometry EL Rank (%)")+
        geom_hline(yintercept =0.5, lty=2, col="grey50")+
        geom_hline(yintercept =2, lty=2, col="grey50")+
        theme(
          # line = element_blank(),
          axis.ticks.length = unit(2, "mm"),
          axis.ticks.y = element_line(color="black"),
          axis.text = element_text(size=14),
          title = element_text(size=12, face = "bold"),
          axis.title = element_text(size=15, face = "bold"))+
        coord_cartesian(ylim = c(0, 11.5))+
        ggtitle(paste0(gene_variant, ":", hla_label))
      
      if(plotting=="on"){ggarrange(a,NULL,b,nrow=1, ncol=3,common.legend = T, legend = "right",align = "hv",
                                   widths = c(1, 0.05, 1)) %>% ggsave(filename=paste0("./netMHCpan/",gene_variant,"_",hla_label,"_netMHCpan.tiff"),
                                                                      device = "tiff", bg="white", width=11,height = 5, units = "in", dpi = 300)}
    }
    
    netMHCpan_rank <- data.frame(ID= netMHCpan$ID, EL_RANK=netMHCpan$EL_Rank, POS = netMHCpan$Pos)
    if(pair ==T & longi ==F){return(b)}
  }
  # mosaic
  table=NULL
  if(DO_LONGI){
    df <- get(paste0("hlaSNPbin_longi_", GENE)) %>% ungroup() %>%
      mutate(
        var = factor(get(paste0(VAR))) %>%  fct_relevel("0","1") %>% fct_recode("present" = "1", "other" = "0") %>% ff_label(paste0("HIV variant ", gene_variant)),
        hla = factor(get(paste0(HLA))) %>%  fct_relevel("0","1") %>% fct_recode("present" = "1","other" = "0") %>% ff_label(paste0("HLA-",hla_label))) %>%
      drop_na(var,hla) %>%
      mutate(var = factor(str_replace_all(var, pattern="present", paste0(variant_label)))%>%  fct_relevel("other"),
             hla = factor(str_replace_all(hla, pattern="present", paste0(hla_label)))%>%  fct_relevel("other")) %>%
      dplyr::select(var,hla, base_uuid, ID, APD,sample.date,time,starts_with("PC"))
    
    df <- merge(df, distinct(all_SHCS %>% ungroup() %>% dplyr::select(ID, SEX)), by="ID") %>% drop_na(var) %>% tibble()
    
    # for labeling the mosaic plot with sample size
    label_numbering = unique(
      df %>% 
        dplyr::mutate(
          hla = hla %>%  fct_relevel("other"),
          var = var %>%  fct_relevel(variant_label)) %>% 
        group_by(.,ID) %>%
        dplyr::arrange(.,(sample.date), .by_group = T) %>% 
        drop_na(.,var,hla) %>% # remove NA
        dplyr::mutate(.,time = as.numeric(difftime(sample.date, min(sample.date), units = "days"))) %>% 
        {if(mosaic_time == "time0") dplyr::filter(.,time==0) else .} %>%
        {if(mosaic_time == "time1")
          filter(.,
                 any(var != variant_label & time == 0) |
                   any(var == variant_label & time != 0),
                 n()>1) %>%
            mutate(.,
                   time0 = min(time),
                   time1 = ifelse(var == variant_label, time, max(time)),
                   time1 = min(time1)) %>%
            filter(.,time == time0 | time == time1) %>% 
            filter(.,time!=0) else .} %>% 
        dplyr::select(var,hla, ID) %>% distinct() %>% with(table(var,hla)))
    
    # mosaic plot with all first sample dates
    plotmosaic <- (df %>% dplyr::mutate(
      hla = hla %>%  fct_relevel("other"),
      var = var %>%  fct_relevel("other",variant_label)) %>% 
        group_by(.,ID) %>%
        dplyr::arrange(.,(sample.date), .by_group = T) %>% 
        drop_na(.,var,hla) %>% # remove NA
        dplyr::mutate(.,time = as.numeric(difftime(sample.date, min(sample.date), units = "days"))) %>% 
        
        {if(mosaic_time == "time0") dplyr::filter(.,time==0) else .} %>%
        {if(mosaic_time == "time1")
          filter(.,
                 any(var != variant_label & time == 0) |
                   any(var == variant_label & time != 0),
                 n()>1) %>%
            mutate(.,
                   time0 = min(time),
                   time1 = ifelse(var == variant_label, time, max(time)),
                   time1 = min(time1)) %>%
            filter(.,time == time0 | time == time1) %>% 
            filter(.,time!=0) else .} %>% 
        dplyr::select(var,hla, ID) %>% distinct() %>% 
        ggplot() +
        geom_mosaic(
          aes(x=product(var,hla),fill=hla),
          divider = ddecker())+
        scale_fill_manual(values = c(colorspace::lighten(color, 0.6),color))+
        annotate("text",label = paste0(label_numbering[1:2]), col = "black", fontface =2, size=4,
                 x = c((label_numbering[1]+label_numbering[2])/(sum(label_numbering))*0.5,
                       (label_numbering[1]+label_numbering[2])/(sum(label_numbering))*0.5),  
                 
                 y = c(1-label_numbering[1]/(label_numbering[1]+label_numbering[2])*0.5,
                       label_numbering[2]/(label_numbering[1]+label_numbering[2])*0.5))+
        
        annotate("text",label = paste0(label_numbering[3:4]), col = "white", fontface =2, size=4,
                 x = c(1-(label_numbering[3]+label_numbering[4])/(sum(label_numbering))*0.5,
                       1-(label_numbering[3]+label_numbering[4])/(sum(label_numbering))*0.5),   
                 
                 y = c(1-label_numbering[3]/(label_numbering[3]+label_numbering[4])*0.5,
                       label_numbering[4]/(label_numbering[3]+label_numbering[4])*0.5))+
        xlab("HLA allele")+
        ylab("HIV AA mutation")+
        theme_blank()+
        theme(
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="None",
          axis.text = element_text(size=12),
          axis.title = element_text(size=13, face = "bold"))) 
    
    if(plotting=="on"){ggsave(plot=plotmosaic,device = "tiff", width = 120, height = 120, bg="white",units="mm", filename = 
                                paste0("./mosaic_longitudinal_",gene_variant,"_",mosaic_time,"_",HLA,".tiff"))}
    
    
    ## Kaplan-Meier survival estimator 
    ## remove those out which already have the mutation and those which have only one sample
    df <- df %>% group_by(ID) %>% 
      dplyr::arrange(.,(sample.date), .by_group = T) %>% 
      drop_na(.,var,hla) %>% # remove NA
      dplyr::mutate(.,time = as.numeric(difftime(sample.date, min(sample.date), units = "days"))) %>% 
      dplyr::filter(.,any(var != variant_label & time == 0)|any(var == variant_label & time != 0)) %>%
      dplyr::filter(.,n()>1) %>% 
      mutate(.,
             time0 = min(time),
             time1 = ifelse(var == variant_label, time, max(time)),
             time1 = min(time1)) %>%
      filter(.,time == time0 | time == time1) %>% 
      filter(time0 != time1)
    
    # check if all have two variables
    df %>% summarise(n= n()) %$% n %>% table()
    
    #remove first sample.date
    df <- df %>% filter(time!=0)
    
    # KM plot
    dependent_os = "Surv(time/365, var)"
    explanatory_os = "hla" # univariable
    # explanatory = c("hla","SEX","PC1","PC2") # multivariable
    
    if(plotting=="on"){tiff(width = 150,height = 120, res=300, units="mm",filename =
                              paste0("./KM_plots/KM_plot_",gene_variant,"_",hla_label,"_confint.tiff"), compression = "none")}
    plot <- df %>% ungroup() %>% dplyr::mutate(var = ifelse(var == "other", FALSE, TRUE),SEX = as.factor(SEX)) %>% dplyr::mutate(hla= relevel(hla,ref="other")) %>% dplyr::select(var,hla,time) %>%
      finalfit::surv_plot(dependent = "Surv(time, var)", explanatory = explanatory_os, 
                          pval = T, pval.coord=c(ceiling(max(df$time)/365-2.5)*365,1.1),
                          pval.method = T, pval.method.coord=c(ceiling(max(df$time)/365-2.5)*365,1.2),
                          xlab="Time (years)", 
                          conf.int = T, conf.int.alpha=0.1,conf.int.style= "ribbon",
                          censor=T,
                          legend.title = "HLA", legend.labs = c("other", hla_label),
                          risk.table = T, 
                          xscale = "d_y", 
                          fun = "cumhaz",
                          ylab=paste0("Cumulative Hazard \nof developing ", gene_variant),
                          ylim=c(0,1.25),
                          break.x.by = (2*365.25),
                          cumevents = F,
                          surv.median.line = "none",
                          cumcensor = F,
                          tables.height = 0.25,risk.table.height	=0.2,
                          palette=c(colorspace::lighten(color, 0.6),color),
                          ncensor.plot=F,
                          linetype=1,
                          ggtheme = theme(axis.ticks.length = unit(2, "mm"), 
                                          panel.background = element_rect(fill="white"),
                                          axis.line = element_line(color="black"),
                                          legend.background = element_rect(fill="white"),
                                          # axis.title.y = element_text(
                                          #   margin = margin(l = 100)),
                                          legend.key = element_rect(fill="white"),
                                          text = element_text(size=12),
                                          legend.title = element_text(size=17),
                                          axis.title.y = element_text(vjust = -3),
                                          legend.text = element_text(size=15),
                                          axis.text = element_text(size=12)))
    
    if(pair ==T & longi ==T){return(plot)}
    print(plot)
    invisible(capture.output(dev.off()))
    
    suppressWarnings(table <- df %>% 
                       ungroup() %>% 
                       mutate(var = ifelse(var == "other", FALSE, TRUE),
                              hla = as.factor(str_remove(pattern = ":",hla)),
                              SEX = as.factor(SEX)) %>%
                       mutate(hla= relevel(hla,ref="other")) %>%
                       finalfit(dependent_os, explanatory_os, add_dependent_label = F) )
  }
  return(list(netMHCpan_rank,table))
}

table_var_change <- full_join(table_var_change,pair_hits %>% dplyr::select(HLA=hla,VAR=var,col)) %>% 
  mutate(Color = ifelse(is.na(col), "black", col)) %>% drop_na(var.n) %>% distinct()

sign_hits_longi <- table_var_change %>% 
  group_by(VAR,HLA,gene) %>%  dplyr::slice(1) %>% 
  ungroup() %>% rowwise() %>% 
  dplyr::mutate(HZ_table = list(longitudinal_analysis(VAR,
                                                      HLA,
                                                      GENE = gene,
                                                      # pair = T,
                                                      plotting = "on",
                                                      DO_LONGI = do_longi,
                                                      # longi = T,
                                                      color=Color,
                                                      mosaic_time = "time1"))) 

sign_hits_longi <- sign_hits_longi %>%unnest_wider(HZ_table, names_sep = ".") %>% dplyr::rename(netMHCpan =HZ_table.1, HZ_table = HZ_table.2)

sign_hits_longi_HZ <- sign_hits_longi %>% dplyr::select(-netMHCpan) %>% unnest(col = c(HZ_table))
sign_hits_longi_HZ <- sign_hits_longi_HZ %>% mutate(
  label = paste0(ifelse(label == "",lag(toupper(label)),toupper(label)))) %>%
  dplyr::rename(uni=`HR (univariable)`,
                multi = `HR (multivariable)`) %>%
  filter(uni != "-") %>% dplyr::select(-all,-levels)

sign_hits_longi_HZ <- sign_hits_longi_HZ %>% 
  dplyr::mutate(
    uni = str_remove_all(uni, pattern="[^[:alnum:][:space:].-]"),
    multi = str_remove_all(multi, pattern="[^[:alnum:][:space:].-]")) %>%
  separate(col = uni, into = c("HZ.uni", "CI.lower.uni","CI.upper.uni", "pval.uni"), sep = " |-") %>%
  separate(col = multi, into = c("HZ.multi", "CI.lower.multi","CI.upper.multi", "pval.multi"), sep = " |-") %>%
  dplyr::mutate(pval.uni = str_remove_all(pval.uni, pattern="[^[:digit:].]")) %>%
  dplyr::mutate(pval.multi = str_remove_all(pval.multi, pattern="[^[:digit:].]"))

sign_hits_longi_HZ <- sign_hits_longi_HZ %>% mutate_at(.vars = colnames(sign_hits_longi_HZ %>% dplyr::select(contains("uni"), contains("multi"))), as.numeric)

sign_hits_Hz <- sign_hits_longi_HZ %>% filter(label == "HLA") %>% drop_na(pval.multi)
sign_hits_Hz <- sign_hits_Hz %>% 
  dplyr::select(!ends_with("uni"), -label)

sign_hits_Hz$gene <- factor(sign_hits_Hz$gene, levels = c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))
sign_hits_longi_netMHCpan <-  sign_hits_longi %>% dplyr::mutate(pos=str_remove_all(VAR, "[^[:digit:]]"))%>%unnest(col = c(netMHCpan)) %>% dplyr::select(-HZ_table) %>% distinct()
sign_hits_longi_netMHCpan$gene <- factor(sign_hits_longi_netMHCpan$gene, levels = c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))
sign_hits_longi_netMHCpan$POS <- as.factor(sign_hits_longi_netMHCpan$POS)
sign_hits_longi_netMHCpan$ID <- as.factor(sign_hits_longi_netMHCpan$ID)
sign_hits_longi_netMHCpan <- sign_hits_longi_netMHCpan %>% group_by(HLA,VAR,gene,POS) %>% arrange(ID,.by_group = T) %>% 
  dplyr::mutate(Rank_log_ratio =log10(EL_RANK / lag(EL_RANK)),
                Rank_ratio =(EL_RANK / lag(EL_RANK)),
                Cons = ifelse(ID=="Consensus", EL_RANK, NA),
                Mut = ifelse(ID=="Mutation", EL_RANK, NA),
  ) %>% 
  dplyr::mutate(Cons = mean(Cons, na.rm=T)) %>% drop_na(Rank_log_ratio) %>% dplyr::select(-ID,-EL_RANK) %>%
  arrange(desc(Rank_log_ratio), pos) %>%
  ungroup()

## thresholds set by NetMHCpan-4.1 (see publication)
sign_hits_longi_netMHCpan <- sign_hits_longi_netMHCpan %>% mutate(
  binding_diff = case_when(
    Cons <= 0.5 & Mut<= 0.5 ~ 'Strong → strong Binder',
    Cons > 0.5 & Cons <= 2 & Mut>0.5 & Mut <= 2 ~ 'Weak → weak Binder',
    
    Cons <= 0.5 & Mut>2 ~ 'Strong → non-Binder',
    Cons > 0.5 & Cons <= 2 & Mut>2 ~ 'Weak → non-Binder',
    
    Cons <= 0.5 & Mut>0.5 & Mut <= 2 ~ 'Strong → weak Binder',
    Cons>0.5 & Cons <= 2 & Mut <= 0.5 ~ 'Weak → strong Binder',
    
    Cons>2 & Mut <= 0.5 ~ 'Non-Binder → strong Binder',
    Cons>2 & Mut>0.5 & Mut <= 2 ~ 'Non-Binder → weak Binder',
    TRUE~NA)
)


sign_hits_Hz <- full_join(sign_hits_Hz, table_var_change) %>% arrange(HLA,VAR,gene)

model_APD_OR <- fread("./hlaSNPbin/model_APD_OR_B.csv") 

#### correlation APD*HLA + longitudinal ####
cor_matrix_APD_longi <- full_join(sign_hits_Hz,model_APD_OR%>% filter(grepl("Interaction", term))) %>% distinct() %>% drop_na(var.n)
cor_matrix_APD_longi$gene <- factor(cor_matrix_APD_longi$gene, levels = c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))
fwrite(cor_matrix_APD_longi, "./cor_matrix_APD_longi.csv")


longi_pair <- pair_hits %>% rowwise() %>% filter(var=="x432R"| var=="x57E") %>% 
  mutate(p1 = list(longitudinal_analysis(VAR=var,
                                         HLA=hla,
                                         GENE = gene,
                                         pair = T,
                                         # longi = T,
                                         DO_LONGI = T,
                                         color=col,
                                         mosaic_time = "time1"))) %>%
  mutate(p2 = list(longitudinal_analysis(VAR=var,
                                         HLA=hla,
                                         GENE = gene,
                                         pair = T,
                                         longi = T,
                                         DO_LONGI = T,
                                         color=col,
                                         mosaic_time = "time1"))) %>% ungroup() 

# Figure 4 ####
arrange_ggsurvplots(longi_pair$p2,axis="r",print = F) %>%
  ggsave(filename=paste0("Figure4BC.tiff"),
         device = "tiff", width=10,height = 4.5, units = "in", dpi = 300)


plot_grid(
  plot_grid(p1, label_size = 15, label_fontfamily = "Arial", label_fontface = "bold",labels = c("A")),
  NULL,
  plot_grid(NULL,NULL, label_size = 15, label_fontfamily = "Arial", label_fontface = "bold",labels = c("B", "C")),
  ncol=1,
  rel_heights = c(1,0.1),
  align = 'hv', axis = 'lr') %>%  
  ggsave(filename=paste0("Figure4ABC.tiff"), 
         device = "tiff",width = 7.3, height = 10, units = "in", dpi = 300)



#### correlation APD*HLA + netMHCpan ####
sign_hits_longi_netMHCpan$binding_diff <- as.factor(sign_hits_longi_netMHCpan$binding_diff)
cor_matrix_APD_netMHCpan <- full_join(sign_hits_longi_netMHCpan, model_APD_OR %>%  filter(term=="Interaction ")%>%dplyr::select(-POS)) %>% drop_na(Rank_log_ratio)
cor_matrix_APD_netMHCpan$binding_diff <- factor(cor_matrix_APD_netMHCpan$binding_diff, levels=c('Non-Binder → weak Binder','Non-Binder → strong Binder',
                                                                                                'Weak → non-Binder',"Weak → weak Binder",'Weak → strong Binder',
                                                                                                'Strong → non-Binder','Strong → weak Binder',"Strong → strong Binder"))
cor_matrix_APD_netMHCpan$gene <- factor(cor_matrix_APD_netMHCpan$gene, levels = c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))
fwrite(cor_matrix_APD_netMHCpan,"./netMHCpan/cor_matrix_APD_netMHCpan.csv")

# Figure 5 ####
p5 <- (cor_matrix_APD_netMHCpan %>% arrange(desc(Color)) %>%
         dplyr::mutate(escape = ifelse(Rank_log_ratio>0 & estimate>1,T,F)) %>%
         dplyr::mutate(gene = factor(firstup(as.character.factor(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) %>%
         
         ggplot(data=., aes(y=Rank_log_ratio, x=log(estimate), color=Color, starshape=gene, alpha=is.na(col)))+
         geom_rect(inherit.aes = F, aes(xmin=0, ymin=0, xmax=Inf, ymax=Inf), alpha=0.2, fill="#FEF3C1")+
         geom_hline(yintercept = 0, color="grey80",linetype=2)+
         geom_vline(xintercept = 0, color="grey80",linetype=2)+
         geom_star(cex=2, fill="white")+
         geom_star(data= . %>% filter(escape),cex=2, aes(fill=Color))+
         scale_starshape_manual(values=c(15, 13, 28, 25, 14, 23, 11)) +
         guides(starshape= guide_legend(label.theme = element_text( ), title="HIV proteins",
                                        override.aes = list(fill = "white"),
                                        title.theme = element_text(face="bold", size = 13)),alpha="none")+
         geom_text_repel(size=3,show.legend = FALSE,max.overlaps = 100,vjust=-1.1,hjust=0.5,
                         aes(label= ifelse(is.na(col), NA, paste0(var,"~",hla, " (",str_remove_all(POS, "[^[:digit:]]"),")"))))+
         scale_color_identity(aesthetics = c("color", "fill"))+
         scale_alpha_manual(values = c(1,0.7))+
         scale_y_continuous(labels = c("-2.0", "-1.0" ,"0.0", "1.0","2.0"))+
         xlab("Log OR APD*HLA")+
         ylab("Impact on binding strength (EL Rank Log Ratio)") +
         annotate("text",size=4,x=Inf, y=2.1, label= "MHC escape mutations", hjust=1.05, vjust=0, color="#BCB375")+
         
         annotate("text",size=4,x=0, y=-Inf, label= "in presence of HLA allele", hjust=0.5, vjust=-1)+
         annotate("text",size=4,x=0.1, y=-Inf, label= "more mutations over time →", hjust=0, vjust=1.7)+
         annotate("text",size=4,x=-0.1, y=-Inf, label= "← fewer mutations", hjust=1, vjust=1.7)+
         
         annotate("text",size=4,x=-Inf, y=0.1, label= "weaker binding upon mutation →", hjust=0,vjust=1.5,angle=90)+
         annotate("text",size=4,x=-Inf, y=-0.1, label= "← stronger binding upon mutation", hjust=1,vjust=1.5,angle=90)+
         
         theme_cowplot()+
         theme(axis.title = element_text(size = 15, face = "bold"),
               axis.text = element_text(size = 12),
               legend.title = element_text(face="bold"),
               axis.ticks.length = unit(2,"mm"))+
         coord_cartesian(ylim=c(-2.2,2), xlim=c(-1.5,1.5),clip="off"))

plot_grid(
  plot_grid(p5, label_size = 20, label_fontfamily = "Arial", label_fontface = "bold",labels = c("A")),
  plot_grid(longi_pair$p1[[1]]+theme(axis.title = element_text(size = 15, face = "bold"),legend.position = "none",
                                     axis.text = element_text(size = 12)),
            longi_pair$p1[[2]]+theme(axis.title = element_text(size = 15, face = "bold"),legend.position = "none",
                                     axis.text = element_text(size = 12)),
            get_legend(longi_pair$p1[[1]]), rel_widths = c(1,1,0.2), ncol=3,
            label_size = 20, label_fontfamily = "Arial", label_fontface = "bold",
            labels = c("B", "C"), align = "hv", axis="l"),
  ncol=1, align = "hv", axis="r",
  rel_heights = c(1,0.65),
  label_size = 20) %>%  
  ggsave(filename=paste0("Figure5.tiff"), bg="white", 
         device = "tiff", width=10,height = 11.5, units = "in", dpi = 300)

# Figure S4 ####
(cor_matrix_APD_netMHCpan %>% 
   mutate(logestimate = log(estimate)) %>% 
   ggplot(data=., aes(x=log(Mut),y=log(Cons),color=logestimate))+
   geom_abline(color="grey65",linetype=1,size=0.7)+
   geom_hline(yintercept=log(0.5), color="grey75", linetype=1)+
   geom_hline(yintercept=log(2), color="grey75", linetype=2)+
   geom_vline(xintercept=log(0.5), color="grey75", linetype=1)+
   geom_vline(xintercept=log(2), color="grey75", linetype=2)+
   
   scale_color_gradientn("HLA*APD OR \ninteraction effect",
                         colours = c("indianred4","red3","indianred3","white","dodgerblue2","blue3","midnightblue"),
                         breaks = log(c(0.25,0.5, 1, 2, 3.5)), ## scale is skewed
                         label= c(0.25,0.5,1.00, 2, 3.5),
                         guide = guide_colorbar(label.position = "right", direction = "vertical",title.position="top",reverse = F,
                                                ticks.colour = "white"))+
   
   geom_point(aes(x=ifelse(binding_diff=="Strong → strong Binder" | binding_diff=="Weak → weak Binder" ,NA,log(Mut)),
                  y=ifelse(binding_diff=="Strong → strong Binder" | binding_diff=="Weak → weak Binder" ,NA,log(Cons))), alpha=1,size=2)+
   
   geom_point(aes(x=ifelse(binding_diff=="Strong → strong Binder" | binding_diff=="Weak → weak Binder" ,log(Mut),NA),
                  y=ifelse(binding_diff=="Strong → strong Binder" | binding_diff=="Weak → weak Binder" ,log(Cons),NA)), alpha=0.3,size=2)+
   
   geom_text_repel(size=2,max.overlaps = 50,vjust=-0.5, aes(label=ifelse(is.na(col),NA, paste0(var,"~",hla, " (",str_remove_all(POS, "[^[:digit:]]"),")"))))+
   
   annotate("text",size=2.5,x=-Inf, y=log(0.5), color="grey10", hjust=1.1,vjust=-1, label= "strong Binder")+
   annotate("text",size=2.5,x=-Inf, y=log(2), color="grey10", hjust=1.1,vjust=-1, label= "weak Binder")+
   annotate("text",size=2.5,x=-Inf, y=Inf, color="grey10", hjust=1.1,vjust=-1, label= "non-Binder")+
   
   annotate("text",size=2.5,y=-Inf, x=log(0.5), color="grey10", hjust=-0.1,vjust=1.1, label= "strong Binder")+
   annotate("text",size=2.5,y=-Inf, x=log(2), color="grey10", hjust=-0.1,vjust=1.1, label= "weak Binder")+
   annotate("text",size=2.5,y=-Inf, x=Inf, color="grey10", hjust=-0.1,vjust=1.1, label= "non-Binder")+
   
   scale_y_reverse(breaks=c(2,0,-2,-4))+
   scale_x_reverse(breaks=c(2,0,-2,-4))+
   xlab("Mutation [log EL Rank]")+
   ylab("Consensus [log EL Rank]")+
   theme_cowplot()+
   theme(legend.box.just = "center",
         axis.title = element_text(size = 12, face = "bold"),
         axis.text = element_text(size = 10),
         legend.title = element_text(face="bold",size=10,hjust = 0),
         # legend.position = "top",
         legend.position = "right",
         
         legend.text = element_text(size=7),
         axis.ticks.length = unit(2,"mm")) +
   coord_cartesian(ylim=c(2.6,-5.2),xlim=c(2.6,-5.2),clip="off"))  %>%
  ggsave(filename="FigureS4.tiff", bg="white", width=7, height=5.5)

