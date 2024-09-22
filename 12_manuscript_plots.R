## manuscript plots from main figures and supplementary figures
## figures from longitudinal or netMHCpan are in 13_longitudinal_netMHC.R file (Figure 4, Figure 5 and Figure S4)
setwd("")
source("packages_load.R")

# Table One ####
tableOne <- meta_combined_db %>% filter(NGS==T) %>% ungroup() %>%
  dplyr::mutate(time_diagn_ART = as.numeric(difftime(ART_START_DATE, date_diagnosis,units = "days"))/365.25)%>% ungroup()

factors= c("SEX", "REGION", "ETHNICITY", "subtype")
variables= c("subtype","ART naive","longi","SEX", "Baseline Age", "REGION", "ETHNICITY","mean log10 RNA", "median CD4", "proviral","standardized APD", "time_diagn_ART")
nonn = c("Baseline Age", "mean log10 RNA", "median CD4","standardized APD", "time_diagn_ART")
filtering=list(
  c("HLA","NGS"),
  c("HLA","NGS","SUBTYPE"),
  c("HLA","NGS","naive"),
  c("HLA","NGS","SUBTYPE","naive"),
  c("HLA","NGS","VL","naive"),
  c("HLA","NGS","SUBTYPE","VL","naive"),
  c("HLA","NGS","longi","naive"),
  c("HLA","NGS","SUBTYPE","longi","naive"))

for(i in 1:length(filtering)){
  filters=filtering[[i]]
  cat("\n",filters)
  print_table <- (tableOne %>% ungroup() %>% drop_na(subtype) %>%
                    # restrict to group of participants needed in analyses A-C
                    filter(NGS == T) %>%
                    {if("HLA" %in% filters) filter(.,HLA == T) else .} %>%
                    {if("VL" %in% filters) filter(.,VL == T) else .} %>%
                    {if("longi" %in% filters) filter(.,longi == T) else .} %>%
                    {if("SUBTYPE" %in% filters) filter(.,subtype == "B") else .} %>%
                    {if("naive" %in% filters) filter(.,`ART naive` ==T) else .} %>%
                    tableone::CreateTableOne(data=., factorVars = factors,vars=variables,addOverall = F,includeNA = T) %>% print(.,nonnormal = nonn))
  write.csv(print_table, paste0("./tableone/TableOne_",str_flatten(filters, collapse = "_"),".csv"))
}

# APD and VL models prep ####
sub="B" # or nonB
proviral = T # TRUE = including proviral in analysis, FALSE = excluding them
sign_hits <- sign_hits_B %>% rbind(., data.frame(
  hla="B_5701",
  var="x242N",
  gene="gag",
  POS=1513,
  FDR=NA,pval=NA))

Fisher <- matrix(nrow=nrow(sign_hits), ncol=7)
colnames(Fisher) <- c("HLA","VAR","gene", "OR","cil","ciu","pval")
OR_APD <- cbind(Fisher[,1:3], coef = NA,interaction_sign=NA,nobs=NA,AIC=NA) %>% as.data.frame()
OR_VL <- as.data.frame(OR_APD) %>% dplyr::rename(r2=AIC)
Fisher <- as.data.frame(Fisher)
ART_naive=T
setwd("")
# loop through pairs
for (i_ in 1:nrow(sign_hits)) {
  vars= sign_hits[i_,]
  position=str_remove_all(vars$NAME, pattern = "[[:letter:]]") 
  reg <- as.character(vars$gene)
  vNAME = vars$var[1]
  vHLA <- vars$hla[1]
  cat("\n",reg,":",vHLA,vNAME,"\n")
  hlaSNPbin <- fread(paste0("./ngs/hlaSNPbin/",ifelse(sub == "nonB", "pooled", sub),"/hlaSNPbin_named_",reg,".csv")) %>% tibble()
  hlaSNPbin <- right_join(NGS_s %>%ungroup()%>% dplyr::select(-fasta,-coverage,-ID ,-ZPHI.ID,-organism,-APD), 
                          hlaSNPbin %>% distinct(), 
                          by="base_uuid") %>% tibble()
  hlaSNPbin$subtype <- fct_infreq(hlaSNPbin$subtype)
  switch <- F
  if(i_ %in% (sign_hits %>%dplyr::mutate(rown= row_number()) %>% group_by(gene) %>% dplyr::slice(1) %$% rown %>% sort())){cat(nrow(hlaSNPbin), "sequences");switch <- T}
  ## filter out proviral seqs
  if (proviral==F){hlaSNPbin<-hlaSNPbin %>% filter(proviral==F)}
  if(switch == T & proviral==F){cat(" minus proviral:", nrow(hlaSNPbin))}
  ## filter out nonB/B seqs
  if (sub == "nonB"){hlaSNPbin<-hlaSNPbin %>% filter(subtype != "B")}
  if(switch == T & sub == "nonB"){cat(" for non-B:", nrow(hlaSNPbin))}
  ## check variant or HLA missing
  df <- tryCatch({hlaSNPbin %>% ungroup() %>% dplyr::select(all_of(vNAME),all_of(vHLA)) %>% drop_na(all_of(vNAME))}, 
                 error = function(e) {message(conditionMessage(e));return(NULL)})
  if(is.null(df)){next}
  # Fisher's exact test
  Fisher$VAR[i_] <-  vNAME
  Fisher$HLA[i_] <- vHLA
  Fisher$gene[i_] <- reg
  Fisher$OR[i_] <- unique(stats::fisher.test(table(df))$estimate)
  Fisher$cil[i_] <- stats::fisher.test(table(df))$conf.int[1]
  Fisher$ciu[i_] <- stats::fisher.test(table(df))$conf.int[2]
  Fisher$pval[i_] <- unique(stats::fisher.test(table(df))$p.value)
  # exclude high APD and sequences under ART - necessary for APD interaction
  hlaSNPbin <- hlaSNPbin %>% filter(APD<0.05) # filter APD<0.05 because of non-linear effect of APD on VL
  hlaSNPbin <- hlaSNPbin %>% filter(ART_naive==T) %>% tibble() # filter if they are ART naive
  if(switch == T ){cat("\nART-naive:", nrow(hlaSNPbin))}
  hlaSNPbin <- hlaSNPbin %>%  ungroup() %>%
    dplyr::mutate(APDz = APD/sd(APD, na.rm = T)) %>% # human PCs (generated on SNP data)
    dplyr::mutate(PC1z = PC1/sd(PC1, na.rm = T),
                  PC2z = PC2/sd(PC2, na.rm = T),
                  PC3z = PC3/sd(PC3, na.rm = T),
                  PC4z = PC4/sd(PC4, na.rm = T),
                  PC5z = PC5/sd(PC5, na.rm = T),
                  PC6z = PC6/sd(PC6, na.rm = T),
                  PC7z = PC7/sd(PC7, na.rm = T),
                  PC8z = PC8/sd(PC8, na.rm = T),
                  PC9z= PC9/sd(PC9, na.rm = T),
                  PC10z = PC10/sd(PC10, na.rm = T)) %>% #viral PCs (Eigenstrat, Nucl whole genome)
    dplyr::mutate(vPC1z = vPC1/sd(vPC1, na.rm = T),
                  vPC2z = vPC2/sd(vPC2, na.rm = T),
                  vPC3z = vPC3/sd(vPC3, na.rm = T),
                  vPC4z = vPC4/sd(vPC4, na.rm = T),
                  vPC5z = vPC5/sd(vPC5, na.rm = T),
                  vPC6z = vPC6/sd(vPC6, na.rm = T),
                  vPC7z = vPC7/sd(vPC7, na.rm = T),
                  vPC8z = vPC8/sd(vPC8, na.rm = T),
                  vPC9z= vPC9/sd(vPC9, na.rm = T),
                  vPC10z = vPC10/sd(vPC10, na.rm = T))
  
  # tests for APD interaction (odds ratios)
  OR_APD$VAR[i_]  <-  vNAME
  OR_APD$HLA[i_] <- vHLA
  OR_APD$gene[i_] <- reg
  hlaSNPbin <- hlaSNPbin %>% dplyr::select(base_uuid,all_of(vNAME), all_of(vHLA), contains("z")) %>% drop_na(all_of(vNAME), all_of(vHLA)) %>%dplyr::mutate_at(.vars = c(vNAME,vHLA), as.factor)
  
  #glm
  glms <- tryCatch({glm(data= hlaSNPbin, formula= paste0(vNAME,"~", as.factor(vHLA), "*APDz+PC1z+PC2z+PC3z+PC4z+PC5z+PC6z+PC7z+PC8z+PC9z+PC10z", ifelse(sub!="B","+vPC1z+vPC2z+vPC3z+vPC4z+vPC5z+vPC6z+vPC7z+vPC8z+vPC9z+vPC10z","")),family = "binomial")}, 
                   warning = function(w) {message("Warning occurred: ", conditionMessage(w));return(NULL)}, error = function(e) {message("An error occurred: ", conditionMessage(e)); return(NULL)})
  if(is.null(glms)){next}
  OR_APD$coef[i_] <-  list(tidy(glms,conf.int = T,conf.level = 0.95))
  OR_APD$interaction_sign[i_] = (OR_APD$coef[[i_]]["p.value"] %>% last()) <=0.05 # filtered only those which have significant interaction and survived Fisher test (HLA~NGS) with FDR<0.2
  OR_APD$nobs[[i_]] = glms %>% nobs()
  OR_APD$AIC[[i_]] <- glms$aic
  # tests for VL (estimates)
  hlaSNPvl <- full_join(fread(paste0("./ngs/phenotypes/SPVL_phenotype_",sub,"_naive.csv")) %>% ungroup() %>% 
                          dplyr::select(-ID,-ZPHI.ID, -APD,-APDz,-starts_with("vPC")),
                        hlaSNPbin, by="base_uuid") %>% tibble() %>% drop_na(mean_logRNA,PC1z)
  # lm
  OR_VL$VAR[i_]  <-  vNAME
  OR_VL$HLA[i_] <- vHLA
  OR_VL$gene[i_] <- reg
  
  lms  <-  tryCatch({lm(data=hlaSNPvl, formula = paste0("mean_logRNA~",vNAME,"*",vHLA,"+APDz+PC1z+PC2z+PC3z+PC4z+PC5z+PC6z+PC7z+PC8z+PC9z+PC10z",  ifelse(sub!="B","+vPC1z+vPC2z+vPC3z+vPC4z+vPC5z+vPC6z+vPC7z+vPC8z+vPC9z+vPC10z","")))}, 
                    warning = function(w) {message("Warning occurred: ", conditionMessage(w));return(NULL)}, error = function(e) {message("An error occurred: ", conditionMessage(e)); return(NULL)})
  if(is.null(lms)){next}
  OR_VL$coef[i_] =  list(tidy(lms,conf.int = T,conf.level = 0.95))
  OR_VL$interaction_sign[i_] = (OR_VL$coef[[i_]]["p.value"] %>% last()) <=0.05 # filtered only those which have significant interaction and survived Fisher test (HLA~NGS) with FDR<0.2
  OR_VL$nobs[[i_]] = lms %>% nobs()
  OR_VL$r2[[i_]] = RsquareAdj(lms)$r.squared
  # changing level of HLA: for APD in presence of HLA
  hlaSNPbin <- hlaSNPbin %>%dplyr::mutate(hla = relevel(get(vHLA),"1"), APDz00 = APDz)
  gl <- glm(data= hlaSNPbin, formula= paste0(vNAME,"~hla*APDz00+PC1z+PC2z+PC3z+PC4z+PC5z+PC6z+PC7z+PC8z+PC9z+PC10z",ifelse(sub!="B","+vPC1z+vPC2z+vPC3z+vPC4z+vPC5z+vPC6z+vPC7z+vPC8z+vPC9z+vPC10z","")),family = "binomial")
  OR_APD$coef[i_] <- list(rbind(data.frame(OR_APD$coef[i_]), tidy(gl,conf.int = T,conf.level = 0.95)[3,]))
}
# safe copy
OR_APD_10 <- OR_APD
OR_VL_10 <- OR_VL

#untangle data
OR_APD <- OR_APD %>% drop_na(interaction_sign) %>% unnest(cols = c(coef,nobs,AIC), names_repair = "unique")%>%group_by(HLA,VAR,gene) #%>% na.omit()
OR_VL <- OR_VL %>% drop_na(interaction_sign) %>% unnest(cols = c(coef,nobs,r2), names_repair = "unique")%>%group_by(HLA,VAR,gene)# %>% na.omit()
Fisher <- tibble(Fisher)%>% na.omit()

## save copy
fwrite(na="NA",as.data.frame(Fisher),file = paste0("./FISHER_",sub,ifelse(proviral,"","_pv"),".csv"), row.names = F)
fwrite(na="NA",as.data.frame(OR_APD),file = paste0("./OR_APD_",sub,ifelse(proviral,"","_pv"),".csv"), row.names = F)
fwrite(na="NA",as.data.frame(OR_VL),file = paste0("./OR_VL_",sub,ifelse(proviral,"","_pv"),".csv"), row.names = F)

## forest plot APD interaction
OR_APD$term <- rep(unlist(c(list("(Intercept)","allele","APD_HLA_absent","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                            ifelse(sub!="B",
                                   list(c("vPC1","vPC2","vPC3","vPC4","vPC5","vPC6","vPC7","vPC8","vPC9","vPC10","Interaction allele:APD","APD_HLA_present")),
                                   list(c("Interaction allele:APD","APD_HLA_present"))))),
                   length(OR_APD$term)/count(OR_APD)[1,4]) 
OR_APD <- OR_APD %>% filter(term !="(Intercept)") %>% filter(!grepl("PC", term)) %>% tibble()
OR_APD <- OR_APD %>% dplyr::select(-std.error, -statistic) %>% group_by(VAR,HLA,gene)

### we have to exp()
OR_APD <- OR_APD %>% 
  dplyr::mutate_at(c("estimate", "conf.low", "conf.high"), exp) %>%
  dplyr::mutate(pval = pval_string(p.value),
                orCI = paste0(sprintf(estimate,fmt = '%#.2f'), " [", sprintf(conf.low ,fmt = '%#.2f'), ", ", sprintf(conf.high ,fmt = '%#.2f'), "]")) %>%
  dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) 


# highest OR top
model_APD_OR <-  OR_APD %>% group_by(VAR,HLA) %>% dplyr::mutate(POS = as.numeric(str_remove_all(VAR, "[[a-z][A-Z]]"))) %>% arrange(gene,POS)

# rename the variables
model_APD_OR <-   model_APD_OR %>% dplyr::mutate(
  hla = insert_colon(str_replace_all(HLA, "_", "*")),
  var = str_remove(VAR, "x"),
  var_names = str_replace(term, "_", " if "),
  var_names = str_replace(var_names, "_", " is "),
  var_names = str_replace(var_names, ":", " x "),
  var_names = str_replace(var_names, "allele", hla),
  var_names = str_remove_all(var_names, "Interaction ")) %>%
  relocate(nobs, .after = interaction_sign)

model_APD_OR <- model_APD_OR %>%dplyr::mutate(term = fct_relevel(term,c("allele", "APD_HLA_absent","APD_HLA_present","Interaction allele:APD"))) %>% arrange(term, .by_group = T)

model_APD_OR$term <- recode(model_APD_OR$term, "allele"= "HLA allele","Interaction allele:APD" = "Interaction ","APD_HLA_absent" = "APD if HLA is absent","APD_HLA_present"= "APD if HLA is present")
model_APD_OR <- model_APD_OR %>% 
  distinct() %>% arrange(gene,POS) 

### VL with NGS*HLA interaction
OR_VL$term <- rep(unlist(c("(Intercept)","variant","allele","APD","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                           ifelse(sub!="B",
                                  list(c("vPC1","vPC2","vPC3","vPC4","vPC5","vPC6","vPC7","vPC8","vPC9","vPC10","Interaction allele:variant")),
                                  list(c("Interaction allele:variant"))))),
                  length(OR_VL$term)/count(OR_VL)[1,4])

OR_VL <- OR_VL %>% filter(term !="(Intercept)") %>% filter(!grepl("PC", term)) %>% dplyr::select(-std.error, -statistic)

OR_VL <- OR_VL %>% rowwise() %>% 
  dplyr::mutate(pval = pval_string(p.value)) %>%
  dplyr::mutate(orCI = paste0(sprintf(estimate,fmt = '%#.2f'), " [", sprintf(conf.low ,fmt = '%#.2f'), ", ", sprintf(conf.high ,fmt = '%#.2f'), "]")) %>% 
  dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) %>%
  group_by(gene, VAR,HLA)

model_VL_B <-  OR_VL %>% dplyr::mutate(POS = as.numeric(str_remove_all(VAR, "[[a-z][A-Z]]"))) %>% arrange(gene,POS)

# rename the variables
model_VL_B <-   model_VL_B %>% dplyr::mutate(
  hla = insert_colon(str_replace_all(HLA, "_", "*")),
  var = str_remove(VAR, "x"),
  var_names = str_replace(term, "_", " if "),
  var_names = str_replace(var_names, "variant", var),
  var_names = str_replace(var_names, ":", " x "),
  var_names = str_replace(var_names, "allele", hla),
  var_names = str_remove_all(var_names, "Interaction "))  %>%
  dplyr::mutate(term=str_remove_all(term, " allele:variant"))

model_VL_B <- model_VL_B %>%dplyr::mutate(term = fct_relevel(term,c("variant", "allele","APD","Interaction variant:allele"))) %>% arrange(term, .by_group = T)
model_VL_B$term <- recode(model_VL_B$term, "variant" = "HIV variant", "allele"= "HLA allele","Interaction variant:allele" = "Interaction")
model_VL_B <- model_VL_B %>% arrange(gene,POS) %>%dplyr::mutate(hlaVar = paste0(HLA,VAR)) %>% filter(hlaVar %in% unique(model_APD_OR %>%dplyr::mutate(hlaVar = paste0(HLA,VAR)) %$% hlaVar)) %>% dplyr::select(-hlaVar)

fwrite(model_APD_OR, paste0("./ngs/hlaSNPbin/model_APD_OR_",ifelse(proviral,"","pv_"),sub,".csv"), na = "NA")
fwrite(model_VL_B, paste0("./ngs/hlaSNPbin/model_VL_B_",ifelse(proviral,"","pv_"),sub,".csv"), na = "NA")

### forest plots ####
sub="B"
model_APD_OR <- fread( paste0("./ngs/hlaSNPbin/model_APD_OR_",ifelse(T,"","pv_"),sub,".csv"), na = "NA")
model_VL_B <- fread( paste0("./ngs/hlaSNPbin/model_VL_B_",ifelse(T,"","pv_"),sub,".csv"), na = "NA")
# model with only the pairs of interest
model_APD_OR2 <- model_APD_OR %>% filter(VAR %in% pair_hits$var)
model_VL_B2 <- model_VL_B %>% filter(VAR %in% pair_hits$var)

source("./forestplot_function.R")
# for (region in unique(sign_hits$gene)) {
region="all"
tryCatch({  
  # for APD-HLA model
  forestploter_table(
    data_name = "model_APD_OR2", 
    batch_size = 9,
    ref_all = T, 
    folder = paste0(" "),
    plotname= "Figure2B.tiff",
    gap=6.5,
    graph.pos=3,
    size=1.3, # not smaller than 0.8
    g_width=1050,
    g_height=1700,
    lineheight = "auto",
    graphwidth = unit(115, 'mm'),
    line.margin=0.16,
    # legend.pos = 1.55,
    show_legend = T,
    gene = region)
  
  # for VL-HLA model
  forestploter_table(
    batch_size = 8,
    data_name = "model_VL_B2", 
    gap=6.5,
    graph.pos=3,
    size=1.3, # not smaller than 0.8
    g_width=1550,
    g_height=1260,
    line.margin=0.09,
    lineheight = "auto",
    show_legend = T,      # if TRUE, the values below have to be adapted
    legend.pos = list(x=1.9,y=0.38),
    graphwidth = unit(115, 'mm'),
    ref_all = T, 
    plotname= "Figure3A.tiff",
    folder = paste0(""),
    gene = region)
  
})
# }

# Figure 2 ####
# ggplot OR Fisher
Fig2 <- list()
Fisher <- Fisher %>% dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) 
ticks <- c(0.05,0.1,0.2,0.5,1,2,5,10,50,250)
l_ticks <- as.character(ticks)

Fig2[[1]] <- (Fisher %>% 
                dplyr::mutate(OR = (OR), cil = (cil),ciu = (ciu)) %>%
                arrange(desc(OR)) %>% drop_na() %>%
                dplyr::mutate(in_pairs_hits = VAR %in% pair_hits$var & HLA %in% pair_hits$hla) %>%
                dplyr::mutate(ID = paste0("ID", row_number())) %>%
                ungroup() %>% rowwise()%>%
                dplyr::mutate(cil = ifelse(cil==0,NA,cil)) %>%
                dplyr::mutate(OR = ifelse(OR==0,NA,OR)) %>%
                dplyr::mutate(ciu = ifelse(ciu==0,NA,ciu)) %>% 
                dplyr::mutate(l_OR = ifelse(cil < min(ticks), min(ticks), NA),
                              u_OR = ifelse(ciu > max(ticks), max(ticks), NA),
                              cil = ifelse(cil < min(ticks), NA, cil),
                              ciu = ifelse(ciu > max(ticks), NA, ciu),
                              NAME= str_remove(VAR, "x"), 
                              hla = insert_colon(str_replace(HLA, "_", "*"))) %>% 
                dplyr::mutate(hlaclass = (str_split(HLA, "_") %>% unlist())[1]) %>%
                
                ggplot(aes(x=(OR),y=reorder(ID, -OR), color=hlaclass,xmax = ciu, xmin = cil,shape= in_pairs_hits))+
                geom_vline(aes(xintercept=1), linewidth = 0.3, linetype="dashed") +
                geom_errorbarh(size = 0.5, height = 1)+
                geom_point() + 
                scale_shape_manual(values = c(15,17)) +
                guides(shape="none",
                       color=guide_legend(override.aes = list(shape=15)))+
                scale_color_manual("HLA gene",
                                   values=c("grey10",'grey40', 'grey65',"red4","red3", 'red1',"indianred2",colorspace::lighten("red",0.5)),
                                   aesthetics = c("color", "fill"))+
                theme_classic()+
                scale_y_discrete(limits = rev) +
                scale_x_continuous(trans="log",
                                   breaks = ticks, 
                                   labels = l_ticks,
                                   limits = c(min(ticks)-min(ticks)/10,max(ticks)))+
                ylab("") +
                theme(axis.title.y=element_blank(),
                      text = element_text(size = 10, face = "bold",),
                      axis.ticks.length.x = unit(2, "mm"),
                      axis.text.x = element_text(size = 10,face = "plain"),
                      axis.title.x = element_text(size = 11, face="bold"),
                      legend.text =  element_text(size = 10,face = "plain"),
                      legend.title =  element_text(size = 11, face="bold"),
                      legend.position = c(x=0.9,y=0.11),
                      axis.text.y=element_blank(), 
                      legend.key.height = unit(5,"mm"),
                      axis.ticks.y =element_blank(), 
                      axis.line.y =element_blank(),
                      strip.placement = "outside")+
                xlab("Odds ratio of the HIV variant to occur")+
                facet_grid(rows = vars(gene), scales = "free_y", space = "free_y", switch="both", shrink = F)+ #facet group
                
                geom_text(hjust=0, x = -Inf,size=2.2,
                          aes(fontface= ifelse(in_pairs_hits, 2, 1),
                              label= ifelse(is.na(OR),NA,paste(str_remove(VAR, "x"), insert_colon(str_replace(HLA, "_", "*")), sep="~"))))+ # add REF
                coord_cartesian(expand = T,clip="off"))

title <- ggdraw() + 
  draw_label(
    "    ",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
Fig2[[1]] <- plot_grid(
  title, Fig2[[1]],
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.03, 1)
)

## APD dynamics on HLA #
pair_hits <- pair_hits %>% dplyr::mutate(gene = factor(firstup(as.character.factor(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) 

Fig2C=list()
#Fig2
filepath= ""

pairs <- pair_hits %>% filter(var=="x432R"|var=="x57E")
for (i in 1:nrow(pairs)) {
  plot=list()
  VAR <- pairs$var[i]
  HLA <- pairs$hla[i]
  variant_label <- str_remove(VAR, "x")
  REGION <- pairs$gene[i]
  
  gene_variant <- paste0(REGION,"-",variant_label)
  hla_label <- insert_colon(str_replace(paste0(HLA), "_", "*"))
  col <-  pairs$col[i]
  hlaSNPbin <- fread(paste0("./ngs/hlaSNPbin/",sub,"/hlaSNPbin_named_",REGION,".csv"))
  mo = (glm(formula= paste0(VAR,"~APD*",HLA),family="binomial",data=hlaSNPbin))
  
  df = hlaSNPbin %>% 
    dplyr::mutate(hla =factor(get(HLA), levels = c("1","0")))%>%
    dplyr::rename(var =all_of(VAR)) %>% tibble() %>%dplyr::select(hla,var,APD) %>% drop_na() %>% dplyr::mutate(p=fitted(mo))
  
  pl <- ggplot(data = df,aes(y=var,x=APD,color= hla,fill=hla))+
    geom_point()+
    stat_smooth(se=T, alpha=0.075, method = "loess",span = 1)+
    stat_smooth(aes(APD, p),se=F,linetype="dashed",size=0.5,method = "glm")+
    xlim(0,0.05)+
    theme_minimal()+
    theme(axis.ticks.length = unit(2, units = "mm"),
          axis.text = element_text(size=10),
          axis.title =  element_text(size=13, face="bold"),
          legend.text = element_text(size=11),
          legend.title = element_text(size=13, face="bold"),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"))+
    scale_y_continuous(breaks=c(0.0,0.5,1.0), labels = c("0.0","0.5","1.0"))+
    guides(fill="none")+
    scale_color_manual(labels=c(hla_label,"other"),
                       values = c(col,colorspace::lighten(col, 0.6)),
                       aesthetics = c("color", "fill")) +
    theme(legend.position ="none")+
    labs(color = "HLA allele",
         x="Average Pairwise Diversity",
         y=paste0("Occurance of ",variant_label))+
    coord_cartesian(ylim = c(0, 1))
  ggsave(pl, filename = paste0(filepath,"APD_", HLA,"_",gene_variant, ".tiff"),width=160, height=125, dpi = 300, units = "mm", device = "tiff", bg="white")
  
  plot[[2]] <- pl
  
  df <- df %>%
    dplyr::mutate(
      var = factor(var) %>%  fct_relevel("0","1") %>% fct_recode("present" = "1", "other" = "0") %>% 
        ff_label(paste0("HIV variant ", gene_variant)),
      hla = factor(hla) %>%  fct_relevel("0","1") %>% fct_recode("present" = "1","other" = "0") %>% 
        ff_label(paste0("HLA-",hla_label))) %>%
    drop_na(var,hla) %>%
    dplyr::mutate(var = factor(str_replace_all(var, pattern="present", paste0(variant_label)))%>%  fct_relevel("other"),
                  hla = factor(str_replace_all(hla, pattern="present", paste0(hla_label)))%>%  fct_relevel("other")) %>%
    dplyr::select(var,hla, APD)
  
  # for labeling the mosaic plot with sample size
  label_numbering = unique(df %>% dplyr::mutate(
    hla = hla %>%  fct_relevel("other"),
    var = var %>%  fct_relevel(variant_label)) %>% with(table(var,hla)))
  
  # mosaic plot with all first sample dates
  pl <- df %>% dplyr::mutate(
    hla = hla %>%  fct_relevel("other"),
    var = var %>%  fct_relevel("other",variant_label)) %>%dplyr::mutate(gene=REGION) %>%
    ggplot() +
    geom_mosaic(
      aes(x=product(var,hla),fill=hla),
      divider = ddecker())+
    scale_fill_manual(values = c(colorspace::lighten(col, 0.6),col))+
    annotate("text",label = paste0(label_numbering[1:2]), col = "black", fontface =2, size=3,
             x = c((label_numbering[1]+label_numbering[2])/(sum(label_numbering))*0.5,
                   (label_numbering[1]+label_numbering[2])/(sum(label_numbering))*0.5),  
             
             y = c(1-label_numbering[1]/(label_numbering[1]+label_numbering[2])*0.5,
                   label_numbering[2]/(label_numbering[1]+label_numbering[2])*0.5))+
    
    annotate("text",label = paste0(label_numbering[3:4]), col = "white", fontface =2, size=3,
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
      axis.text = element_text(size=11, color = "grey40"),
      axis.title = element_text(size=12, face = "bold", color="black"))
  
  ggsave(pl, device = "tiff", width = 100, height = 95, bg="white",units="mm", 
         filename = paste0(filepath,"mosaic_", HLA,"_",gene_variant, ".tiff"))
  plot[[1]] <- pl
  # align mosaic and scatterplot
  title <- ggdraw() + 
    draw_label(
      paste0(firstup(REGION)),
      fontface = 'bold',
      x = 0,
      hjust = 0
    )  +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 50)
    )
  
  arranged <- (plot_grid(title,plot_grid(plotlist =  plot,axis = "tblr"
                                         # ,rel_widths = c(1.7, 2)
  ),
  ncol=1, rel_heights = c(0.09,1))) 
  Fig2C[[i]] <- arranged
  
}
Fig2[[3]] <- plot_grid(plotlist =Fig2C, ncol=1)

## scatterplot APD OR 
df_scatterplot <- full_join(model_APD_OR, pair_hits %>%dplyr::select(HLA=hla, VAR=var,gene,col)) %>% filter(grepl(term, pattern="APD"))
df_scatterplot <- df_scatterplot %>%
  dplyr::select(HLA,VAR,gene,interaction_sign,POS,term,estimate,col,nobs,hla,var) %>% 
  pivot_wider(names_from = "term", values_from = estimate)%>% 
  dplyr::mutate(`Variant Pair` = paste0(gene," - ",VAR," ~ ", HLA, " (n=",nobs,")")) %>%
  dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) %>%
  ungroup()%>%
  dplyr::mutate(estimate_y= ifelse(`APD if HLA is present`>summary(`APD if HLA is present`)[5]+1.5, as.numeric(summary(`APD if HLA is present`)[5])+1.5, `APD if HLA is present`)) %>%
  dplyr::mutate(estimate_x= ifelse(`APD if HLA is absent`>summary(`APD if HLA is absent`)[5]+1.5, as.numeric(summary(`APD if HLA is absent`)[5])+1.5, `APD if HLA is absent`)) %>%
  dplyr::mutate(Color = ifelse(is.na(col), "black", col)) 


Fig2[[4]] <- (df_scatterplot %>%
                ggplot(aes(x=estimate_x,y=estimate_y, 
                           color=Color, fill=Color,
                           starshape=gene,
                           alpha= is.na(col)))+
                
                scale_x_continuous(trans = "log",breaks=c(0.5,1,2))+
                scale_y_continuous(trans = "log",breaks=c(0.3,0.5,1,2,3.5))+
                
                geom_rect(xmin =-Inf, xmax=0,ymin=0,ymax=Inf, alpha=0.2, fill="#FEF3C1",color="#FEF3C1")+
                
                geom_text(data=.%>%filter(estimate_y==max(estimate_y)),color="grey40", alpha=0.9,
                          aes(label="*", hjust=-1))+
                
                geom_hline(yintercept = 1, col="grey10",linetype="dashed", size=0.4)+
                geom_vline(xintercept = 1,col="grey10",linetype="dashed", size=0.4) +
                
                coord_cartesian(clip="off",ylim = c(0.3,3.5),xlim = c(exp(-0.9),exp(0.9)))+
                
                scale_color_identity(aesthetics = c("fill","color"))+
                guides(color="none",fill="none",alpha="none") + 
                scale_alpha_discrete(range = c(1, 0.6))+
                labs(x="OR of APD if allele is absent",y="OR of APD if allele is present")+
                
                geom_star(size=2.5, fill="white")+
                geom_star(data = . %>% filter(!is.na(col)),size=3)+
                
                scale_starshape_manual(values=c(15, 13, 28, 25, 12, 14, 23, 11))+
                guides(starshape= guide_legend(label.theme = element_text(face="plain",size = 10), title="HIV protein",
                                               override.aes = list(fill = "grey96",size=2),
                                               ncol=2,byrow = T,
                                               title.theme = element_text(face="bold", size = 12))) +
                
                geom_text_repel(size=4,max.overlaps = 10, 
                                vjust=1.1,hjust=-0.5,
                                aes(label= ifelse(is.na(col),NA,paste(var,hla, sep="~"))))+
                theme_classic()+
                theme(axis.ticks.length = unit(2, "mm"),
                      axis.ticks = element_line(linewidth = 0.3), 
                      panel.grid = element_blank(),
                      text = element_text(size=13),
                      legend.text = element_text(size=13),
                      legend.position = c(0.12,0.14),
                      legend.title = element_text(size=14),
                      legend.key.size = unit(2.5,"mm"),
                      axis.title = element_text(size=14, face = "bold"))) 

title <- ggdraw() + 
  draw_label(
    paste0("     "),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  draw_label(
    paste0("  "),
    fontface = 'bold.italic',
    x = 0,
    hjust = -1.85
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 50)
  )

Fig2[[4]] <- plot_grid(title,plot_grid(NULL,Fig2[[4]],rel_widths = c(0.05,1)),ncol=1,rel_heights = c(0.05,1))

plot_grid(plotlist =Fig2, ncol=2, rel_heights = c(2,1.1),rel_widths = c(1,1.1), axis="t",align = "hv",labels = "AUTO", label_size = 16) %>%
  ggsave(filename = "./Figure2ACD.tiff",
         width = 30, height = 40, units = "cm", dpi = 300,device = "tiff")

# Figure 3 ####
## boxplot violinplot 
filepath = " "
plots <- VL_violin_boxplot(loop_df = unique(pair_hits) %>% filter(var=="x57E"|var=="x432R"), filepath = " ",
                           pair="all", w=120, h=100)

title <- ggdraw() + 
  draw_label(
    paste0("sample size"),
    x = -1,y=0.955,
    size=11,
    hjust = 0) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 50))


### scatterplot VL
# model_VL_B from above
VL_scatterplot <- model_VL_B  %>% 
  ungroup() %>% 
  dplyr::select(-conf.low,- conf.high,-nobs,-orCI, -HLA,-VAR, -var_names,-pval,-r2,-POS) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  dplyr::select(-APD) %>%
  group_by(hla,var,gene) %>%
  dplyr::summarise(gene,interaction_sign,hla,var,`HIV variant`=mean(`HIV variant`,na.rm=T), 
                   `HLA allele`=mean(`HLA allele`,na.rm=T),Interaction=mean(Interaction,na.rm=T)) %>% distinct()

VL_scatterplot$interaction_sign <- factor(VL_scatterplot$interaction_sign, levels = c(F,T))
VL_scatterplot <- full_join(VL_scatterplot,pair_hits %>%
                              dplyr::mutate(var= str_remove_all(pattern = "x", var),
                                            hla=insert_colon(str_replace_all(pattern="_", replacement = "*", hla))) %>%
                              dplyr::select(var,gene,hla,col)) %>% 
  dplyr::mutate(Color = ifelse(is.na(col), "black", col))  %>%
  arrange(desc(Color),interaction_sign, desc(Interaction))

plots[[3]] <-  
  (VL_scatterplot %>% drop_na(interaction_sign) %>% ungroup() %>% 
     dplyr::mutate(col=ifelse(var=="9L", "black", col)) %>%
     dplyr::mutate(gene = factor(firstup(as.character.factor(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) %>%
     
     ggplot(aes(x = `HIV variant`,y=Interaction,color=Color,fill=interaction_sign,starshape=gene,alpha= is.na(col))) +
     geom_hline(yintercept = 0, col="black", linetype="dashed")+
     geom_vline(xintercept = 0,col="black", linetype="dashed")+
     scale_color_identity(aesthetics = c("color"))+
     theme_classic()+ # HIV AA variant | HLA allele | Interaction
     ylab("HIV:HLA interaction effect on VL")+
     xlab("HIV AA variant")+
     scale_alpha_discrete(range = c(1, 0.7))+
     scale_fill_manual("Interaction", values = c("white", "black"),labels = c("n.s.", "p<0.05"))+
     scale_x_continuous(n.breaks = 3)+
     geom_star(size=2.1)+
     scale_starshape_manual(values=c(15, 13, 28, 25, 12, 14, 23, 11))+
     guides(fill = guide_legend(reverse = TRUE,
                                label.theme = element_text(face="plain"),
                                override.aes = list(fill = c("black","white"), starshape=c(15,15)),
                                title.theme = element_text(face="bold", size = 13)),
            alpha = "none",
            starshape= guide_legend(label.theme = element_text(face="italic"), title="HIV genes",
                                    override.aes = list(fill = "white"),
                                    # ncol=2,
                                    title.theme = element_text(face="bold", size = 13))) +
     geom_star(data=. %>% filter(!is.na(col)) %>% filter(interaction_sign==T),
               size=2.2,
               fill=VL_scatterplot %>% 
                 dplyr::mutate(col=ifelse(var=="9L", "black", col)) %>%
                 
                 filter(!is.na(col)) %>% 
                 filter(interaction_sign==T) %$% col)+
     
     theme(axis.ticks.length = unit(2, "mm"),
           axis.text = element_text(size=12),
           legend.position = "none",
           axis.title = element_text(size=13, face = "bold"))+
     geom_text_repel(size=3.5,max.overlaps = 10, vjust=-2,hjust=0.5,
                     aes(label= ifelse(
                       is.na(col)|
                         interaction_sign==F,
                       NA,paste(var,hla, sep="~"))))+
     coord_cartesian(xlim = c(-0.5,0.5), ylim=c(-0.5,1.1)))


plots[[4]] <- 
  (VL_scatterplot %>% drop_na(interaction_sign) %>% ungroup() %>% 
     dplyr::mutate(col=ifelse(var=="9L", "black", col)) %>%
     dplyr::mutate(gene = factor(firstup(as.character.factor(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) %>%
     ggplot(aes(x = `HLA allele`,y=Interaction, color=Color,fill=interaction_sign,starshape=gene,alpha= is.na(col))) +
     geom_rect(xmin = -Inf, xmax=0,ymin=0, ymax=Inf, alpha=0.2, fill="#FEF3C1",color="#FEF3C1")+
     geom_hline(yintercept = 0, col="black", linetype="dashed")+
     geom_vline(xintercept = 0,col="black", linetype="dashed")+
     scale_color_identity(aesthetics = c("color"))+
     theme_classic()+ # HIV AA variant | HLA allele | Interaction
     xlab("HLA allele")+
     ylab("HIV:HLA interaction effect on VL")+
     scale_alpha_discrete(range = c(1, 0.7))+
     scale_fill_manual("Interaction", values = c("white", "black"),labels = c("n.s.", "p<0.05"))+
     # scale_shape_manual("HIV genes",values=c(21,22,23,7,9,12,24,25))+
     theme(axis.ticks.length = unit(2, "mm"),
           axis.text = element_text(size=12),
           axis.title = element_text(size=13, face = "bold"))+
     #special symbols
     geom_star(size=2.1)+
     scale_starshape_manual(values=c(15, 13, 28, 25, 12, 14, 23, 11))+
     guides(fill = guide_legend(reverse = TRUE,
                                label.theme = element_text(face="plain"),
                                override.aes = list(fill = c("black","white"), starshape=c(15,15)),
                                # ncol=2,
                                title.vjust=1,
                                title.theme = element_text(face="bold", size = 13)),
            alpha = "none",
            starshape= guide_legend(
              label.theme = element_text(), 
              title="HIV proteins",
              title.vjust=1,
              override.aes = list(fill = "white"),
              # nrow=2,byrow = T,
              title.theme = element_text(face="bold", size = 13))) +
     geom_star(data=. %>% filter(!is.na(col)) %>% filter(interaction_sign==T),
               size=2.2,
               fill=VL_scatterplot %>% 
                 dplyr::mutate(col=ifelse(var=="9L", "black", col)) %>%
                 filter(!is.na(col)) %>% filter(interaction_sign==T) %$% col)+
     geom_text_repel(size=3.5,max.overlaps = 10, vjust=-1.1,hjust=1,
                     aes(label= ifelse(
                       is.na(col)|
                         interaction_sign==F,
                       NA,paste(var,hla, sep="~"))))+
     coord_cartesian(xlim = c(-1.1,0.45), ylim=c(-0.5,1.1)))

legend <- get_legend(plots[[4]] +
                       theme(legend.position = "right",
                             legend.title.align = 0,
                             legend.direction = "vertical"))

plots[[4]] <- plots[[4]] + theme(legend.position = "none")

plot_grid(
  plot_grid(plots[[1]], plots[[2]],title,rel_widths = c(1,1,0.3), align = "h",rel_heights = c(1,1,0.08),nrow=1,label_size = 16, labels = c("B", "C")),
  plot_grid(plots[[3]], plots[[4]],legend,label_size = 16, labels = c("D", "E"), align = "v", axis = 'l',rel_widths = c(1,1,0.3), nrow=1),
  ncol=1,
  align = "h", axis = "b")%>%
  ggsave(filename = paste0(filepath,"Figure3BCDE.tiff"),
         width=300, height = 210, units="mm", dpi = 300,device = "tiff")

# Figure S2 ####
# complete df with all OR of Fisher test and GLM 
## can chose between interaction (*) or simple multivariate (+) model
rmlike("Fisher_comb_")
neg = pos = tot = count=0
sub="B"
region_names= c("vif", "vpr", "tat","rev","vpu","nef", "gag", "env", "pol")
for (region in region_names){
  cat("\n",region, " ")
  
  hlaSNPbin <- fread(paste0("./ngs/hlaSNPbin/",sub,"/hlaSNPbin_named_",region,".csv")) %>% tibble()
  
  hlaSNPbin <- left_join(hlaSNPbin %>% distinct(), NGS_s %>%ungroup()%>% dplyr::select(base_uuid,subtype,ART_naive, starts_with("vPC")),by = join_by(base_uuid)) %>% tibble()
  hlaSNPbin$subtype <- relevel(as.factor(hlaSNPbin$subtype), ref=ifelse(sub=="B", sub, as.factor(hlaSNPbin$subtype) %>% table() %>% sort() %>% last() %>% names()))
  
  # necessary for *APD
  hlaSNPbin <- hlaSNPbin %>% filter(APD<0.05) # filter APD<0.05 because of non-linear effect of APD on VL
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
    dplyr::mutate(vPC1z = vPC1/sd(vPC1, na.rm = T),
                  vPC2z = vPC2/sd(vPC2, na.rm = T),
                  vPC3z = vPC3/sd(vPC3, na.rm = T),
                  vPC4z = vPC4/sd(vPC4, na.rm = T),
                  vPC5z = vPC5/sd(vPC5, na.rm = T),
                  vPC6z = vPC6/sd(vPC6, na.rm = T),
                  vPC7z = vPC7/sd(vPC7, na.rm = T),
                  vPC8z = vPC8/sd(vPC8, na.rm = T),
                  vPC9z= vPC9/sd(vPC9, na.rm = T),
                  vPC10z = vPC10/sd(vPC10, na.rm = T))
  
  translation <- fread(paste0("./ngs/hlaSNPbin/",sub,"/translation_named_", region,".csv")) %>% tibble()
  Fisher <- fread(paste0("./ngs/hlaSNPbin/B/Fisher_named_",region,".csv"))
  Fisher_comb <- Fisher %>% pivot_longer(cols = !hla_allele, names_to = "var", values_to = "FDR") %>% drop_na(FDR) %>% distinct()
  
  count = count+nrow(Fisher_comb %>% filter(!var %in% names(hlaSNPbin)))
  Fisher_comb <- Fisher_comb %>% filter(var %in% names(hlaSNPbin))
  tot=tot+as.numeric(nrow(Fisher_comb)) # count the Fisher FDR<0.2
  
  Fisher_comb <- Fisher_comb %>% rowwise() %>% 
    dplyr::mutate(glms = ## add +APD for the descriptive part of 459
                    list(glm(data= hlaSNPbin, formula= (paste0(var,"~",hla_allele,
                                                               # "*APDz",
                                                               "+PC1z+PC2z+PC3z+PC4z+PC5z+PC6z+PC7z+PC8z+PC9z+PC10z")),
                             family = "binomial"))) %>%
    dplyr::mutate(coef = list(tidy(glms)),
                  # ci_APD = list(confint(glms)), # this step needs a long time, so if not needed -> comment
                  nobs_glm = list(glms %>% nobs())) %>%  
    unnest(cols = c(coef,
                    # ci_APD,
                    nobs_glm), names_repair = "unique") %>% 
    dplyr::select(-glms)
  
  neg = neg+Fisher_comb %>% filter(estimate<0 & p.value<0.05 & grepl("_", term)  & !grepl(":", term)) %>% nrow() %>% as.numeric()
  pos= pos+Fisher_comb %>% filter(estimate>0 & p.value<0.05 & grepl("_", term)  & !grepl(":", term)) %>% nrow()%>% as.numeric()
  
  Fisher_comb <- Fisher_comb %>% 
    filter(p.value<=0.05 & 
             grepl("_", term)  &
             !grepl(":", term)) # for multivariate (+) model
  # grepl(":", term)) # filter only those which have significant association and survived Fisher test (HLA~NGS) with FDR<0.2
  
  hlaSNPbin[grep(x = colnames(hlaSNPbin), pattern = "x|_")] <- hlaSNPbin[grep(x = colnames(hlaSNPbin), pattern = "x|_")] %>%dplyr::mutate_if(is.numeric, as.factor)
  if(nrow(Fisher_comb) == 0){cat("none", fill=T); next}
  
  for (k in 1:nrow(Fisher_comb)) {
    df <- hlaSNPbin %>% ungroup() %>% 
      dplyr::select(var = Fisher_comb$var[k], hla = Fisher_comb$hla_allele[k])
    Fisher_comb$OR[k] <- unique(stats::fisher.test(table(df))$estimate)
    Fisher_comb$ci_lower[k] <- stats::fisher.test(table(df))$conf.int[1]
    Fisher_comb$ci_uper[k] <- stats::fisher.test(table(df))$conf.int[2]
    Fisher_comb$pval[k] <- unique(stats::fisher.test(table(df))$p.value)
    
    Fisher_comb$id[k] = paste0(as.data.frame(table(df),responseName="Freq") %>% arrange(desc(Freq)) %>% 
                                 dplyr::select(-Freq)  %>% dplyr::mutate(combi = paste0(var,hla))  %>% 
                                 dplyr::mutate(combi = paste(combi, collapse = "")) %>% dplyr::select(combi)%>% dplyr::slice(1))}
  
  Fisher_comb <- Fisher_comb %>% arrange(OR)
  
  Fisher_comb <- Fisher_comb %>% group_by(id) %>% arrange(id) %>% dplyr::mutate(gene = region) %>% ungroup()
  assign(paste0("Fisher_comb_", region), Fisher_comb)
}

cat(tot,"total\n", pos, "positive assoc\n", neg, "negative assoc \n", count, "missed bc different dataset" )
rm(Fisher_comb)
Fisher_comb <- do.call(rbind, mget(ls()[grepl("Fisher_comb_*", ls())])) #%>% filter(grepl(term, pattern=":"))
Fisher_comb$gene <- factor(Fisher_comb$gene, levels = c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))

# Apply p.adjust() function to the filtered columns
Fisher_comb <- Fisher_comb %>% dplyr::mutate(pval.adj = p.adjust (pval, method='fdr', n=nrow(.)))
Fisher_comb <- Fisher_comb %>% dplyr::select(-FDR,-term,-estimate,-std.error,-statistic,-p.value)
Fisher_comb <- Fisher_comb %>% dplyr::mutate(gene = factor(firstup(as.character.factor(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef")))) 

ticks <- c(0.01,0.03,0.1,0.2,0.5,1,2,5,10,30,100)
l_ticks <- as.character(ticks)

(Fisher_comb %>% ungroup() %>% tibble() %>% filter(pval.adj<0.05) %>%
    arrange(desc(OR)) %>% 
    dplyr::mutate(OR = (OR), ci_lower = (ci_lower),ci_uper = (ci_uper)) %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="A") %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="A") %>% ## for adding space in the plot, not changing any data
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="A") %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="B") %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="A") %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="A") %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="A") %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="B") %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="A") %>%
    add_row(.after = "Vpu", gene=factor("Vpu"), hla_allele="A") %>%
    dplyr::mutate(ID = paste0("ID", row_number())) %>%
    ungroup() %>% rowwise()%>%
    dplyr::mutate(HLA_class = (str_split(hla_allele, "_") %>% unlist())[1]) %>%
    dplyr::mutate(l_OR = ifelse(ci_lower < min(ticks), min(ticks), NA),
                  l_ci = ifelse(ci_lower < min(ticks), NA, ci_lower),
                  NAME= str_remove(var, "x"), 
                  HLA = str_replace(hla_allele, "_", "*")) %>% 
    
    ggplot(aes(x=OR,y=reorder(ID, -OR), color=HLA_class,xmax = ci_uper, xmin = l_ci))+
    coord_cartesian(xlim = c(min(ticks), max(Fisher_comb$ci_uper)),expand = T,clip="off")+
    geom_errorbarh(size = 0.5, height = 2)+
    geom_segment(aes(x= OR,xend = l_OR, yend=reorder(ID, -OR)), position = position_nudge(x = 0.005, y = 0),arrow = arrow(length = unit(1, "mm")),show.legend = F)+
    geom_errorbarh(size = 0.5, height = 2, aes(xmax = ci_uper, xmin = OR))+
    geom_point(size = 1) + theme_classic()+
    geom_vline(aes(xintercept=1), size = 0.5, linetype="dashed") +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(trans="log",
                       breaks = ticks,
                       labels = l_ticks) +
    
    guides(color=guide_legend(title = "HLA gene"))+
    ylab("") +
    xlab("Odds ratio of the HIV variant to occur")+
    facet_grid(rows = vars(gene), scales = "free_y", space = "free" , switch="both", shrink = F)+ #facet group
    scale_color_manual(values=c("grey10",'grey50', 'grey80',"red4","red3", 'red1',"indianred2",colorspace::lighten("red",0.5)))+
    theme(axis.title.y=element_blank(),
          text = element_text(size = 13.5),
          axis.ticks.length.x = unit(2, "mm"),
          axis.text.x = element_text(size = 13,face = "plain"),
          axis.title.x = element_text(size = 14, face="bold"),
          legend.text =  element_text(size = 17,face = "italic"),
          legend.title =  element_text(size = 18, face="bold"),
          legend.position = c(x=0.9,y=0.13),
          axis.text.y=element_blank(), 
          axis.ticks.y =element_blank(), 
          axis.line.y =element_blank(),
          strip.placement = "outside")) %>%
  ggsave(filename = "./FigureS2.tiff", 
         width = 25, height = 30, units = "cm", bg="white", dpi = 300)

# Figure S3 ####
# HLA-wide association analysis on VL 
# merge HLA and VL data
SPVL_hlaBin <- full_join(SPVL_30 %>%dplyr::select(ID, mean_logRNA,SEX), hlaBin, by="ID") %>% drop_na(PC1, mean_logRNA)  %>% distinct()

# generate linear regression models for all alleles
df <- data.frame(HLA = names(hlaBin %>%dplyr::select(contains("_"),-rs333_T))) %>% rowwise() %>%
  dplyr::mutate(HLA_class = (str_split(HLA, "_") %>% unlist())[1])

df <- df %>% rowwise() %>%
  dplyr::mutate(lms = list(lm(data= SPVL_hlaBin,formula=paste0("mean_logRNA~",HLA,"+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")))) %>%
  dplyr::mutate(coef =  list(tidy(lms)),
                nobs = list(lms %>% nobs())) %>%  
  unnest(cols = c(coef,nobs), names_repair = "unique") %>% 
  dplyr::select(-lms) %>% filter(grepl(term, pattern="_")) %>%dplyr::select(HLA, HLA_class, estimate, p.value)

df <- df %>%dplyr::mutate(log10_pval = -log10(p.value)) %>% rowwise() %>%
  dplyr::mutate(HLA = as.factor(HLA)) %>%
  dplyr::mutate(label_hla = as.factor(insert_colon((str_split(HLA, "_") %>% unlist())[2])),
                HLA_class = as.factor (HLA_class),
                HLA_class_wrong = str_remove_all(HLA, pattern="[^[:alpha:]]"),
                hla_2 = paste0(HLA_class_wrong,"*",(str_split(label_hla, ":") %>% unlist())[1]),
                hla_1 = paste0(str_split(label_hla, ":") %>% unlist())[1])

(df %>% ungroup() %>%dplyr::mutate(POS= row_number()) %>% arrange(POS) %>%    
    dplyr::mutate(sign= case_when(p.value> 0.05 ~ "0.5",
                                  p.value< 0.05/nrow(df) ~ "1",
                                  TRUE ~ "0.75")) %>%
    ungroup() %>% rowwise()%>% 
    ggplot(aes(y=log10_pval,x=hla_1))+
    geom_point(aes(alpha=sign))+
    geom_hline(show.legend = T, aes(yintercept = -log10(0.05),linetype="p=0.05"))+
    geom_hline(show.legend = T, aes(yintercept = -log10(0.05/nrow(df)),linetype="p=0.05/n"))+
    theme_classic()+
    scale_color_manual("HLA gene",values= c('#800000',"#332268",'#4363d8','#911eb4', '#f58249', '#008080',
                                            '#3cb44b', "#EFB939","#FF6600","#CC0066", "#731941"))+
    scale_alpha_manual(values = c(0.4,0.7,1))+
    scale_linetype_manual("Thresholds",values=c(2,3))+
    coord_cartesian(ylim= c(0.3,13))+
    facet_grid(cols = vars(HLA_class), scales = "free", space = "free", switch="both", shrink = T)+ #facet group
    guides(alpha="none",
           linetype=guide_legend(reverse = T),
           color= guide_legend(override.aes = aes(label = "", linetype = 0)))+
    geom_text_repel(size=3,vjust=1.05,
                    aes(label=ifelse(log10_pval> -log10(0.05/nrow(df)), paste0(HLA_class,"*",label_hla), NA)))+
    labs(x= "HLA alleles",
         y= expression(bold(~-log ["10"]~ "p-value")))+
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.5, size=6,face ="plain"),
          axis.ticks.length = unit(1, "mm"),
          axis.text = element_text(size=10),
          legend.title =  element_text(size=13,face ="bold"),
          axis.title = element_text(size=13,face ="bold"),
          strip.text = element_text(size = 7, face = "italic"),
          panel.spacing = unit(-0.08, "lines"),
          strip.placement = "outside")) %>%
  ggsave(width = 27, height = 11,bg="white",dpi = 300,
         units = "cm",filename = "./FigureS3.tiff", device = "tiff")

# Figure S5 ####
#### correlation APD VL #### 
FigS5 <- list()

B_OR <- fread(paste0("./ngs/hlaSNPbin/model_APD_OR_B.csv")) %>% tibble()%>% filter(var!="242N") %>%
  dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))))

B_VL <- fread(paste0("./ngs/hlaSNPbin/model_VL_B_B.csv")) %>% tibble()%>% filter(var!="242N") %>%
  dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))))

FigS5[[1]] <- (full_join(B_VL %>% filter(term=="Interaction") %>% dplyr::select(HLA, VAR, gene,estimate_VL=estimate,interaction_sign), 
                         B_OR %>% filter(term=="Interaction ") %>% dplyr::select(HLA, VAR, gene,estimate_APD=estimate),
                         by = join_by(HLA, VAR, gene)) %>% filter(VAR!="x242N") %>% 
                 ggplot(aes(y=estimate_VL,x=estimate_APD,color=interaction_sign))+
                 geom_hline(yintercept = 0, linetype="dashed")+
                 geom_vline(xintercept = 1, linetype="dashed")+
                 geom_point()+
                 geom_point(data=.%>%filter(interaction_sign==T), color="black")+
                 scale_color_manual("Interaction on VL", values = c("grey70", "black"),labels = c("all", "p<0.05"),aesthetics = c("color", "linetype"))+
                 theme_classic()+
                 labs(y="Impact of HLA*HIV interaction on VL",
                      x="Impact of HLA*APD interaction on HIV")+
                 scale_x_continuous(trans = "log", breaks=c(0.1,0.25,0.5,1,2,4),label=c(0.1,0.25,0.5,1,2,4), limits=c(0.2,5.1))+
                 scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1),label=c(-1,-0.5,0,0.5,1), limits=c(-1.2,1.2))+
                 stat_cor(method="spearman",cor.coef.name = "rho", data=(.%>%filter(interaction_sign==T)), color="black", 
                          hjust=0, vjust=0,r.digits = 2,p.accuracy =0.001)+
                 stat_cor(method="spearman",cor.coef.name = "rho",color="grey70", 
                          hjust=0, vjust=1,r.digits = 2,p.accuracy = 0.001)+
                 stat_smooth(method="lm", color="grey70", alpha=0.07,fullrange = T,size=0.5)+
                 stat_smooth(method="lm", data=.%>%filter(interaction_sign==T), color="black", alpha=0.08,fullrange = T,size=0.5)+
                 
                 guides(color = guide_legend(reverse = T, size=3,
                                             label.theme = element_text(size = 11),
                                             title.theme = element_text(face="bold", size = 13))) +
                 
                 theme(axis.title = element_text(face="bold", size = 14),
                       axis.text = element_text(size=12),
                       legend.position = c(0.88,0.08),
                       legend.text.align = 0,
                       legend.key.height = unit(4,"mm"),
                       legend.key.size = unit(10,"mm"),
                       axis.ticks.length = unit(2,"mm"))) 

full_join(B_VL %>% filter(term=="Interaction") %>% dplyr::select(HLA, VAR, gene,estimate_VL=estimate,interaction_sign), 
          B_OR %>% filter(term=="Interaction ") %>% dplyr::select(HLA, VAR, gene,estimate_APD=estimate),
          by = join_by(HLA, VAR, gene)) %>% filter(VAR!="x242N") %>% ungroup() %>% 
  with(cor.test(y=estimate_VL,x=estimate_APD,method = "spearman"))   

#### correlation APD HR ####
# from 02_longitudinal_netMHCpan
cor_matrix_APD_longi <- fread("./longitudinal/cor_matrix_APD_longi.csv")

FigS5[[2]] <-  (cor_matrix_APD_longi %>% filter(VAR!="x242N") %>% tibble()  %>%
                  dplyr::mutate(HZ_sign = pval.multi <= 0.05) %>% 
                  dplyr::mutate(HZ.multi = ifelse(is.infinite(CI.upper.multi)|CI.lower.multi==0,NA,HZ.multi)) %>%
                  ggplot(aes(y=HZ.multi,x=estimate))+
                  geom_hline(yintercept = 1, linetype="dashed")+
                  geom_vline(xintercept = 1, linetype="dashed")+
                  geom_point(aes(color=HZ_sign))+
                  geom_point(data=.%>%filter(HZ_sign==T), color="black")+
                  scale_color_manual("Hazard Ratio", values = c("grey70", "black"),labels = c("n.s.", "p<0.05"),aesthetics = c("color"))+
                  theme_classic()+
                  labs(y="HIV variant Occurence over time (Hazard Ratio)",
                       x="Impact of HLA*APD interaction on HIV")+
                  scale_x_continuous(trans = "log", breaks=c(0.1,0.25,0.5,1,2,4),label=c(0.1,0.25,0.5,1,2,4), limits=c(0.2,5.1))+
                  scale_y_continuous(trans = "log", breaks=c(1e-1,0.25,0.5,1,2,4,1e1,25),label=c(1e-1,0.25,0.5,1,2,4,1e1,25))+
                  guides(color = guide_legend(reverse = T, size=3,
                                              label.theme = element_text(size = 11),
                                              title.theme = element_text(face="bold", size = 13))) +
                  stat_cor(method="spearman", color="grey30",cor.coef.name = "rho", hjust=0, vjust=0,r.digits = 2,p.accuracy =0.001,aes(y=HZ.multi,x=estimate))+
                  stat_smooth(method="lm", color="grey30", alpha=0.07,fullrange = T,size=0.5,aes(y=HZ.multi,x=estimate),inherit.aes = F)+
                  
                  theme(axis.title = element_text(face="bold", size = 14),
                        axis.text = element_text(size=12),
                        legend.position = c(0.88,0.08),
                        legend.text.align = 0,
                        legend.key.height = unit(4,"mm"),
                        legend.key.size = unit(10,"mm"),
                        axis.ticks.length = unit(2,"mm"))) 

cor_matrix_APD_longi %>% filter(VAR!="x242N") %>%      
  dplyr::mutate(HZ_sign = pval.multi <= 0.05) %>% 
  dplyr::mutate(HZ.multi = ifelse(is.infinite(CI.upper.multi)|CI.lower.multi==0,NA,HZ.multi)) %>%
  with(cor.test(y=HZ.multi,x=estimate,method = "spearman"))

#### correlation APD ELrank ####
# from 02_longitudinal_netMHCpan
cor_matrix_APD_netMHCpan <- fread("./netMHCpan/cor_matrix_APD_netMHCpan.csv")

FigS5[[3]] <-  (cor_matrix_APD_netMHCpan %>% filter(VAR!="x242N") %>%     
                  dplyr::mutate(escape_pos = ifelse(Rank_log_ratio>0 & estimate>1,T,F)) %>% 
                  dplyr::mutate(escape_neg = ifelse(Rank_log_ratio<0 & estimate<1,T,F)) %>% 
                  dplyr::mutate(escape = (escape_pos|escape_neg) & binding_diff!="Strong → strong Binder"&binding_diff!="Weak → weak Binder") %>%
                  group_by(HLA,VAR,gene) %>%
                  arrange(desc(abs(Rank_log_ratio)), .by_group = T) %>% dplyr::slice(1) %>% ungroup() %>%
                  ggplot(aes(y=Rank_ratio,x=estimate))+
                  geom_hline(yintercept = 1, linetype="dashed")+
                  geom_vline(xintercept = 1, linetype="dashed")+
                  geom_point(aes(color=escape))+
                  geom_point(data=.%>%filter(escape==T), color="black")+
                  scale_color_manual("NetMHCpan EL Rank", values = c("grey70", "black"),labels = c("all", "change in category"),aesthetics = c("color"))+
                  theme_classic()+
                  labs(y="EL Rank Ratio (Mutation/Consensus)",
                       x="Impact of HLA*APD interaction on HIV")+
                  scale_x_continuous(trans = "log", breaks=c(0.1,0.25,0.5,1,2,4),label=c(0.1,0.25,0.5,1,2,4), limits=c(0.2,5.1))+
                  scale_y_continuous(trans = "log", breaks=c(1e-2,1e-1,0.5,1,2,1e1,1e2,1e3),labels = c(1e-2,1e-1,0.5,1,2,1e1,1e2,1e3), limits = c(1e-2,1e2))+
                  
                  guides(color = guide_legend(reverse = T, size=3,
                                              label.theme = element_text(size = 11),
                                              title.theme = element_text(face="bold", size = 13))) +
                  
                  stat_cor(method="spearman",cor.coef.name = "rho", color="black", data=.%>%filter(escape==T), inherit.aes = F, aes(y=Rank_ratio,x=estimate),
                           hjust=0, vjust=0,r.digits = 2,p.accuracy =0.001)+
                  stat_cor(method="spearman",cor.coef.name = "rho", color="grey70", inherit.aes = F, aes(y=Rank_ratio,x=estimate),
                           hjust=0, vjust=1,r.digits = 2,p.accuracy =0.001)+
                  stat_smooth(method="lm", color="black", alpha=0.07,fullrange = T,size=0.5,data=.%>%filter(escape==T))+
                  stat_smooth(method="lm", color="grey70", alpha=0.07,fullrange = T,size=0.5)+
                  theme(axis.title = element_text(face="bold", size = 14),
                        axis.text = element_text(size=12),
                        legend.position = c(0.85,0.08),
                        legend.text.align = 0,
                        legend.key.height = unit(4,"mm"),
                        legend.key.size = unit(10,"mm"),
                        axis.ticks.length = unit(2,"mm"))) 

# direction of the effects in comparison to APD*HLA direction
tibble(cor_matrix_APD_netMHCpan) %>% filter(VAR!="x242N") %>%     
  dplyr::mutate(escape_pos = ifelse(Rank_log_ratio>0 & estimate>1,T,F)) %>% 
  dplyr::mutate(escape_neg = ifelse(Rank_log_ratio<0 & estimate<1,T,F)) %>%
  dplyr::mutate(escape = (escape_pos|escape_neg) & binding_diff!="Strong → strong Binder"&binding_diff!="Weak → weak Binder") %>%
  group_by(HLA,VAR,gene) %>%
  arrange(desc(abs(Rank_log_ratio)), .by_group = T) %>% 
  # filter(escape) %>%
  # dplyr::slice(1)  %>%
  filter(!any(escape_neg)) %>%
  filter(!any(escape_pos))

## change in binding category
tibble(cor_matrix_APD_netMHCpan) %>% filter(VAR!="x242N") %>% 
  dplyr::mutate(escape_pos = ifelse(Rank_log_ratio>0 & estimate>1,T,F)) %>% 
  dplyr::mutate(escape_neg = ifelse(Rank_log_ratio<0 & estimate<1,T,F)) %>%
  dplyr::mutate(escape = (escape_pos|escape_neg) & binding_diff!="Strong → strong Binder"&binding_diff!="Weak → weak Binder") %>%
  group_by(HLA,VAR,gene) %>%
  arrange(desc(abs(Rank_log_ratio)), .by_group = T) %>% 
  ### epitopes
  # with(paste("stronger upon mutation:", sum(binding_diff=="Weak → strong Binder"|binding_diff=="Non-Binder → weak Binder"|binding_diff=="Non-Binder → strong Binder"), "epitopes"))
  # with(paste("weaker upon mutation:", sum(binding_diff=="Strong → weak Binder"|binding_diff=="Weak → non-Binder"|binding_diff=="Strong → non-Binder"), "epitopes")) 
  # with(paste("same upon mutation:", sum(binding_diff=="Strong → strong Binder"|binding_diff=="Weak → weak Binder"), "epitopes"))
  
  ### pairs
  # filter(binding_diff=="Weak → strong Binder"|binding_diff=="Non-Binder → weak Binder"|binding_diff=="Non-Binder → strong Binder") %>% count() %$% paste(nrow(.),"pairs with epitopes stronger upon mutation")
  # filter(binding_diff=="Strong → weak Binder"|binding_diff=="Weak → non-Binder"|binding_diff=="Strong → non-Binder") %>% count() %$% paste(nrow(.),"pairs with epitopes weaker upon mutation")
  filter(binding_diff=="Strong → strong Binder"|binding_diff=="Weak → weak Binder") %>% count() %$% paste(nrow(.),"pairs with epitopes that don't change category upon mutation")

#### correlation B vs non-B #######
nonB_OR <- fread(paste0("./ngs/hlaSNPbin/model_APD_OR_nonB",".csv")) %>% tibble()%>% filter(var!="242N") %>%
  dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))))

FigS5[[5]] <- (full_join(B_OR,nonB_OR, by=c("HLA","VAR","gene","term","var_names","POS","hla","var"),suffix = c(".B",".nonB")) %>%group_by(HLA,VAR,gene) %>% 
                 dplyr::select(-contains("AIC"),-contains("X2.5"),-contains("nobs"),-contains("X97.5"),-contains("ci")) %>% 
                 drop_na() %>% filter(term == "Interaction ") %>% arrange(estimate.nonB) %>%
                 ggplot(., aes(x=estimate.B, y=estimate.nonB,color=interaction_sign.nonB))+
                 geom_hline(yintercept = 1, linetype="dashed")+
                 geom_vline(xintercept = 1, linetype="dashed")+
                 geom_point()+
                 geom_point(data = . %>% filter(interaction_sign.nonB))+
                 
                 scale_x_continuous(trans = "log", breaks=c(0.1,0.25,0.5,1,2,4),label=c(0.1,0.25,0.5,1,2,4), limits=c(0.2,5.1))+
                 scale_y_continuous(trans = "log", breaks=c(1e-2,1e-1,0.5,1,2,1e1,1e2,1e3),labels = c(1e-2,1e-1,0.5,1,2,1e1,1e2,1e3), limits = c(1e-2,1e2))+
                 stat_cor(vjust=0, color="black",method="spearman",cor.coef.name = "rho",r.digits = 2,p.accuracy =0.001)+
                 stat_smooth(method="lm", color="black", alpha=0.08, size=0.5,fullrange = T)+
                 guides(color= guide_legend(override.aes = aes(label = "", linetype = 0),
                                            reverse = T, size=3,
                                            label.theme = element_text(size = 11),
                                            title.theme = element_text(face="bold", size = 13)))+
                 
                 scale_color_manual(aesthetics = c("fill", "color"), "Interaction (non-B)", labels=c("n.s.", "p<0.05"),values=c("grey70", "black"))+
                 theme_classic()+
                 xlab("Impact of HLA*APD interaction on HIV\nSubtype B")+
                 ylab("Subtype non-B\nImpact of HLA*APD interaction on HIV")+
                 theme(axis.title = element_text(face="bold", size = 14),
                       axis.text = element_text(size=12),
                       legend.position = c(0.85,0.08),
                       legend.text.align = 0,
                       legend.key.height = unit(4,"mm"),
                       legend.key.size = unit(10,"mm"),
                       axis.ticks.length = unit(2,"mm"))+
                 
                 geom_text_repel(size=2.5,hjust=1.1,vjust=0.5,
                                 aes(label= ifelse(interaction_sign.nonB,paste0(gene,var,"~",hla),NA)))  # add REF
) 

pooled_OR <- fread(paste0("./ngs/hlaSNPbin/model_APD_OR_pooled",".csv")) %>% tibble()%>% filter(var!="242N") %>%
  dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))))

# pooled #
FigS5[[6]] <- (full_join(B_OR,pooled_OR, by=c("HLA","VAR","gene","term","var_names","POS","hla","var"),suffix = c(".B",".pooled")) %>%group_by(HLA,VAR,gene) %>% 
                 dplyr::select(-contains("AIC"),-contains("X2.5"),-contains("nobs"),-contains("X97.5"),-contains("ci")) %>% 
                 drop_na() %>%
                 filter(term == "Interaction ") %>% arrange(estimate.pooled) %>% 
                 ggplot(., aes(x=estimate.B, y=estimate.pooled,color=interaction_sign.pooled))+
                 geom_hline(yintercept = 1, linetype="dashed")+
                 geom_vline(xintercept = 1, linetype="dashed")+
                 geom_point()+
                 geom_point(data = . %>% filter(interaction_sign.pooled))+
                 
                 scale_x_continuous(trans = "log", breaks=c(0.1,0.25,0.5,1,2,4),label=c(0.1,0.25,0.5,1,2,4), limits=c(0.2,5.1))+
                 scale_y_continuous(trans = "log", breaks=c(1e-1,0.2,0.5,1,2,5,1e1),labels = c(1e-1,0.2,0.5,1,2,5,1e1), limits = c(1e-1,1e1))+
                 stat_cor(vjust=0, color="black",method="spearman",cor.coef.name = "rho",r.digits = 2,p.accuracy =0.001)+
                 stat_smooth(method="lm", color="black", alpha=0.08, size=0.5,fullrange = T)+
                 guides(color= guide_legend(override.aes = aes(label = "", linetype = 0),
                                            reverse = T, size=3,
                                            label.theme = element_text(size = 11),
                                            title.theme = element_text(face="bold", size = 13)))+
                 
                 scale_color_manual(aesthetics = c("fill", "color"), "Interaction (non-B & B)", labels=c("p<0.05","n.s."),values=c("black","grey70"))+
                 theme_classic()+
                 xlab("Impact of HLA*APD interaction on HIV\nSubtype B")+
                 ylab("Subtype non-B and B pooled\nImpact of HLA*APD interaction on HIV")+
                 theme(axis.title = element_text(face="bold", size = 14),
                       axis.text = element_text(size=12),
                       legend.title = element_text(face="bold", size = 13),
                       legend.text = element_text(size = 11),
                       legend.position = c(0.85,0.08),
                       legend.text.align = 0,
                       legend.key.height = unit(4,"mm"),
                       legend.key.size = unit(10,"mm"),
                       axis.ticks.length = unit(2,"mm"))
                 
) 

#### correlation proviral excluded vs included #######
B_npv_OR <- fread(paste0("./ngs/hlaSNPbin/model_APD_OR_",ifelse(F,"","pv_"),"B",".csv")) %>% tibble()%>% filter(var!="242N")%>%
  dplyr::mutate(gene = factor(firstup(as.character(gene)),levels=firstup(c("gag","pol","vif","vpr","tat","rev","vpu","env","nef"))))

FigS5[[4]] <- (full_join(B_OR,B_npv_OR, by=c("HLA","VAR","gene","term","var_names","POS","hla","var"),suffix = c(".pv",".nonpv")) %>%group_by(HLA,VAR,gene) %>% 
                 dplyr::select(-contains("AIC"),-contains("X2.5"),-contains("nobs"),-contains("X97.5"),-contains("ci")) %>% 
                 drop_na() %>% filter(term == "Interaction ") %>% arrange(estimate.nonpv) %>%
                 dplyr::mutate(ratio = estimate.pv/estimate.nonpv) %>%
                 ggplot(., 
                        aes(x=(estimate.pv), y=(estimate.nonpv),
                            color=interaction_sign.nonpv))+
                 geom_hline(yintercept = 1, linetype="dashed")+
                 geom_vline(xintercept = 1, linetype="dashed")+
                 geom_point()+
                 geom_point(data = . %>% filter(interaction_sign.nonpv))+
                 scale_x_continuous(trans = "log", breaks=c(0.1,0.25,0.5,1,2,4),label=c(0.1,0.25,0.5,1,2,4), limits=c(0.2,5.1))+
                 scale_y_continuous(trans = "log", breaks=c(0.1,0.25,0.5,1,2,4),label=c(0.1,0.25,0.5,1,2,4), limits=c(0.2,5))+
                 stat_cor(vjust=0, color="black",method="spearman",cor.coef.name = "rho",r.digits = 2,p.accuracy =0.001)+
                 stat_smooth(method="lm", color="black", alpha=0.08, size=0.5,fullrange = T)+
                 guides(color= guide_legend(override.aes = aes(label = "", linetype = 0)))+
                 
                 scale_color_manual(aesthetics = c("fill", "color"), "Interaction (excl.)", labels=c("n.s.", "p<0.05"),values=c("grey70", "black"))+
                 theme_classic()+
                 guides(color= guide_legend(override.aes = aes(label = "", linetype = 0),
                                            reverse = T, size=3,
                                            label.theme = element_text(size = 11),
                                            title.theme = element_text(face="bold", size = 13)))+
                 xlab("Impact of HLA*APD interaction on HIV\nProviral included")+
                 ylab("Proviral excluded\nImpact of HLA*APD interaction on HIV")+
                 theme(axis.title = element_text(face="bold", size = 14),
                       axis.text = element_text(size=12),
                       legend.position = c(0.87,0.08),
                       legend.text.align = 0,
                       legend.key.height = unit(4,"mm"),
                       legend.key.size = unit(10,"mm"),
                       axis.ticks.length = unit(2,"mm"))+
                 
                 geom_text_repel(size=2.3,vjust=1.1,hjust=1.5,point.padding = 1,max.overlaps = 20,
                                 aes(label= ifelse(interaction_sign.nonpv==F,paste0(gene,var,"~",hla),NA)))  # add REF
) 

plot_grid(plotlist = FigS5, labels = "AUTO",label_size = 16, ncol=2, align="hv") %>%
  ggsave(filename = paste0("~/switchdrive/qid/masterThesisNNJ/manuscript figures/FigureS5.tiff"), 
         width = 15, height = 18, bg="white",dpi = 300) ## 26 is the optimal width for word

#### ORs from weaker vs stronger binding upon mutation ####
weak <- sign_hits_longi_netMHCpan %>% 
  group_by(HLA,VAR,gene) %>%
  # arrange(desc(abs(Rank_log_ratio)), .by_group = T) %>% dplyr::slice(1) %>% ungroup() %>%
  dplyr::filter(binding_diff=="Strong → non-Binder"|binding_diff=="Weak → non-Binder"|binding_diff=="Strong → weak Binder") %>%dplyr::mutate(hlaVar=paste0(HLA,VAR))
strong <- sign_hits_longi_netMHCpan %>%  group_by(HLA,VAR,gene) %>%
  # arrange(desc(abs(Rank_log_ratio)), .by_group = T) %>% dplyr::slice(1) %>% ungroup() %>%
  filter(binding_diff=="Non-Binder → weak Binder"|binding_diff=="Weak → strong Binder"|binding_diff=="Non-Binder → strong Binder")%>%dplyr::mutate(hlaVar=paste0(HLA,VAR))

model_APD_OR %>% group_by(HLA,VAR,gene) %>% filter(term=="Interaction " & VAR!="x242N") %>%
  dplyr::mutate(hlaVar=paste0(HLA,VAR)) %>%
  dplyr::mutate(binding_change = factor(case_when(hlaVar %in% weak$hlaVar~ "weaker",
                                                  hlaVar %in% strong$hlaVar ~ "stronger",
                                                  TRUE ~ NA))) %>% with(table(binding_change,exclude=F))
drop_na(binding_change) %>%
  ggplot(aes(y=(estimate), x=binding_change))+
  geom_boxplot() + xlab("binding upon mutation")+ylab("Interaction APD*HLA on viral variant (OR)")+
  geom_point(position = position_jitter(width=0.15))+
  stat_compare_means(method = "t.test")

model_APD_OR %>% group_by(HLA,VAR,gene) %>% filter(term=="Interaction " & VAR!="x242N") %>%
  dplyr::mutate(hlaVar=paste0(HLA,VAR)) %>%
  dplyr::mutate(binding_change = factor(case_when(hlaVar %in% weak$hlaVar~ "weaker",
                                                  hlaVar %in% strong$hlaVar ~ "stronger",
                                                  TRUE ~ NA))) %>% drop_na(binding_change) %>% 
  # group_by(binding_change) %>% count()
  with(t.test((estimate)~binding_change))


# Figure S6 ####
#### linkage  
## between HLA on same viral variants
linkage_sign_hits <- model_APD_OR %>% filter(term=="Interaction ") %>% group_by(VAR) %>% filter(n()>1) %>%dplyr::mutate(region=(gene))
j=0
summary_genes <- data.frame(linkage_sign_hits)%>% dplyr::select(VAR, region) %>% distinct() %>%dplyr::mutate(model=NA)

for (reg in unique(summary_genes$region)) {
  hlaSNPbin <- fread(paste0("./ngs/hlaSNPbin/B/hlaSNPbin_named_",reg,".csv")) %>% tibble()
  
  hlaSNPbin <- right_join(NGS_s %>%ungroup()%>% dplyr::select(-fasta,-coverage,-ID ,-ZPHI.ID,-organism,-APD), 
                          hlaSNPbin %>% distinct(), 
                          by="base_uuid") %>% tibble()
  hlaSNPbin$subtype <- fct_infreq(hlaSNPbin$subtype)
  
  # necessary for *APD
  hlaSNPbin <- hlaSNPbin %>% filter(APD<0.05) # filter APD<0.05 because of non-linear effect of APD on VL
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
    dplyr::mutate(vPC1z = vPC1/sd(vPC1, na.rm = T),
                  vPC2z = vPC2/sd(vPC2, na.rm = T),
                  vPC3z = vPC3/sd(vPC3, na.rm = T),
                  vPC4z = vPC4/sd(vPC4, na.rm = T),
                  vPC5z = vPC5/sd(vPC5, na.rm = T),
                  vPC6z = vPC6/sd(vPC6, na.rm = T),
                  vPC7z = vPC7/sd(vPC7, na.rm = T),
                  vPC8z = vPC8/sd(vPC8, na.rm = T),
                  vPC9z= vPC9/sd(vPC9, na.rm = T),
                  vPC10z = vPC10/sd(vPC10, na.rm = T))
  
  hlaSNPbin[grep(x = colnames(hlaSNPbin), pattern = "x|_")] <- hlaSNPbin[grep(x = colnames(hlaSNPbin), pattern = "x|_")] %>%dplyr::mutate_if(is.numeric, as.factor)
  
  hits <- linkage_sign_hits %>% filter(region==reg)  %>% group_by(VAR)
  
  m=list()
  for (i in 1:nrow(hits %>% count)) {
    variant = unique(hits %$%VAR)[i]
    hla1= (hits %>% filter(VAR==variant) %$% HLA) [1]
    hla2= (hits %>% filter(VAR==variant) %$% HLA) [2]
    hla3= (hits %>% filter(VAR==variant) %$% HLA) [3]
    hla4 = (hits %>% filter(VAR==variant) %$% HLA) [4]
    
    m <-glm(data= hlaSNPbin, formula= paste0(variant,"~", 
                                             hla1, "*APDz+", 
                                             hla2, "*APDz+", 
                                             ifelse(is.na(hla3), "", paste0(hla3, "*APDz+")), 
                                             ifelse(is.na(hla4), "", paste0(hla4, "*APDz+")), 
                                             "PC1z+PC2z+PC3z+PC4z+PC5z+PC6z+PC7z+PC8z+PC9z+PC10z"),
            family = "binomial") %>% tidy(conf.int = T,conf.level = 0.95) %>% filter(!grepl(term, pattern="PC|Intercept"))
    j=j+1
    summary_genes[j,]$model  <-  list(m)
  }
}

summary_model <- summary_genes %>% unnest(model) %>%dplyr::select(-std.error,-statistic) %>% filter(grepl(term, pattern=":")) 

summary_model <- summary_model %>%dplyr::mutate(HLA= str_remove_all(term, pattern="APDz|:")) %>%dplyr::mutate(term= paste0(HLA, ":APDz"))
summary_model <- left_join(
  summary_model,
  model_APD_OR %>% ungroup() %>% filter(term=="Interaction ") %>%dplyr::mutate(term= paste0(HLA,1, ":APDz"),HLA=paste0(HLA,1), region = gene),
  by = c("VAR", "HLA", "term", "region")) %>% tibble() 
color_palette <- c(
  "black",'grey0',"grey5","grey10",
  "grey15","grey20","grey25", "grey30",
  "grey35","grey40","grey45",
  "grey50","grey55","grey60","grey65",
  "grey70","grey75","grey80","cornsilk3", "red4","red3",  'red1',"indianred2")

summary_model <- summary_model %>%dplyr::mutate(VAR_region=paste0(firstup(region), str_remove_all(VAR,"x"))) %>%dplyr::mutate(shapes="linkage") %>%arrange(region,POS)

summary_model$shapes <- factor(summary_model$shapes, levels = c("linkage", "original"))
(summary_model %>%
    dplyr::mutate(hla=as.factor(hla),
                  VAR_region=factor(VAR_region,levels=rev(summary_model %$% VAR_region %>% unique))) %>%
    add_row(VAR_region=factor("Gag26R"),POS=11, region=factor("vpr"), shapes="original") %>%
    dplyr::mutate_at(.vars = c("estimate.x", "conf.low.x", "conf.high.x"), exp) %>%
    
    ggplot(aes(x=(estimate.x),color=hla, y=VAR_region, shape=shapes))+
    geom_errorbarh(aes(xmin=(conf.low.x), xmax=(conf.high.x)),size = 0.5, height = 0.5,position=position_dodge(width = 0.5))+
    geom_point(size=2.5,position=position_dodge(width = 0.5))+
    geom_point(data=.%>%dplyr::mutate(estimate.x=estimate.y)%>%mutate(shapes="original"), 
               size=1.5, shape=24, position=position_dodge(width = 0.5)) + # original data added
    geom_text(aes(x=0.037,color=hla, y=VAR_region,label=hla),inherit.aes = F, size=3,hjust=0,position=position_dodge(width = 0.6))+
    scale_shape_manual("Model",values=c(16,24),label=c("with linkage", "without linkage"))+
    geom_vline(xintercept = 1,linetype=2) +
    scale_x_continuous(trans = "log", breaks=c(5e-2,1e-1,0.5,1,2,1e1,1e2,5e1),labels = c(5e-2,1e-1,0.5,1,2,1e1,1e2,5e1))+
    scale_color_manual("HLA alleles",
                       values=color_palette)+
    coord_cartesian(clip="off", xlim = c(5e-2,5e1))+
    theme_cowplot()+
    guides(color="none")+
    theme(axis.text.y = element_text(hjust=0),
          axis.title = element_text(face="bold",size=16),
          axis.ticks.length = unit(2, "mm"),
          legend.position = c(0.85, 0.1),
          axis.title.y = element_text(vjust=1.5))+
    labs(x="Impact of HLA*APD interaction on HIV variants",
         y="HIV variants")) %>%
  ggsave(.,filename="./FigureS6.tiff", bg="white", height = 9.5, width = 12)

#### Table S1 ####
# table of all 98 pairs - Supplementary Table 1
TabS1 <- full_join(
  model_APD_OR %>% filter(term=="Interaction ")%>%dplyr::mutate(`HIV variant ~ HLA allele` = paste0(firstup(gene),var,"~",hla)) %>% dplyr::select(HLA, VAR,gene,`HIV variant ~ HLA allele`,estimate_APD=orCI),
  model_VL_B %>% filter(term=="Interaction") %>% dplyr::select(HLA, VAR, gene,estimate_VL=orCI,interaction_sign), 
  by = join_by(HLA, VAR, gene)) %>% filter(VAR!="x242N") 

TabS1 <- full_join(TabS1,sign_hits_Hz %>%dplyr::mutate(sign_pval=pval.multi<0.05,
                                                       estimate_HR=paste0(sprintf(HZ.multi,fmt = '%#.2f'), " [", sprintf(CI.lower.multi ,fmt = '%#.2f'), ", ", sprintf(CI.upper.multi ,fmt = '%#.2f'), "]"))%>%
                     dplyr::select(HLA,VAR,gene,estimate_HR,sign_pval))%>% filter(VAR!="x242N") 

TabS1 <- full_join(TabS1,
                   sign_hits_longi_netMHCpan %>% 
                     dplyr::mutate(change=binding_diff!="Strong → strong Binder"&binding_diff!="Weak → weak Binder",
                                   estimate_MHC=paste0("pos",str_remove_all(POS,"[[:alpha:]]" ),": ",Cons, " → ",Mut)) %>% 
                     group_by(HLA,VAR,gene) %>%
                     dplyr::mutate(estimate_MHC=paste(estimate_MHC,collapse = "\n"),
                                   change=paste(change,collapse = "\n")) %>% ungroup() %>%
                     dplyr::select(HLA,VAR,gene,estimate_MHC,change))%>% filter(VAR!="x242N") %>% distinct()

fwrite(TabS1 %>% ungroup() %>% 
         dplyr::rename(
           `Interaction APD*HLA on HIV variant\nOR [95% CI]` = estimate_APD,
           `Interaction HIV*HLA on VL\nβ [95% CI]` = estimate_VL,
           `Time-to-event presence/absence HLA\nHR [95% CI]` = estimate_HR,
           `NetMHCpan EL Rank \nPosition: Consensus → Mutation` = estimate_MHC) %>%
         dplyr::select(-HLA,-VAR,-gene),
       file = "./TableS1.csv")

## published pairs ### for table S1
published_pairs <- read_excel("journal.pone.0006687.s001.xls") ## https://www.hiv.lanl.gov/content/immunology/pdf/2010/escape_article_supplement.html
published_pairs <- published_pairs %>%dplyr::mutate(Protein = ifelse(Protein=="INT", "IN", Protein))%>%dplyr::mutate(gene_variant= paste0(Protein, Codon, `Amino Acid`)) 
interactionAPD_model_output <- fread("./interactionAPD_model_output.csv")[,1:9] %>% tibble() %>%
  rowwise() %>%
  dplyr::mutate(hla_short = str_remove_all(HLA, "[^[:alpha:][:digit:]]")) %>%
  dplyr::mutate(hla_short = str_sub(hla_short,end = 3L))

published_pairs  %>% filter(gene_variant %in% toupper(interactionAPD_model_output$VAR_region) | gene_variant %in% toupper(interactionAPD_model_output$VAR_gene)) %>% 
  filter(HLA %in% interactionAPD_model_output$hla_short) %>% arrange(`Direct or Indirect`) %>%
  dplyr::mutate(gene_variant = paste0(gene_variant,"~",HLA)) %>% View()

interactionAPD_model_output %>%dplyr::mutate(published = ((VAR_region %in% published_pairs$gene_variant | VAR_gene %in% published_pairs$gene_variant) & hla_short %in% published_pairs$HLA)) %>% 
  filter(DIRECTION=="direct") %>% View()
