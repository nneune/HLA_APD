# packages

require(pacman)
require(tidyverse)
pacman::p_load(ape,readstata13,ggstar,
               stats,vegan,foreign,lme4,lmerTest,gdata,readxl,Hmisc,ggplot2,ggforce,survival,stringr,survminer,ggfortify,readr,lubridate,generics,viridis,hablar,
               purrr,reshape2,data.table,svMisc,sjPlot,sjlabelled,sjmisc,cowplot,forestplot,ggmosaic,RColorBrewer,pixiedust,xfun,finalfit,doParallel,foreach,
               Biostrings,DECIPHER,
               forestploter,ggrepel,tableone)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")


# load forestplot function
source("./forestplot_function.R")

"rmlike" <- function(...) {
  names <- sapply(
    match.call(expand.dots = FALSE)$..., as.character)
  names = paste(names,collapse="|")
  Vars <- ls(1)
  r <- Vars[grep(paste("^(",names,").*",sep=""),Vars)]
  rm(list=r,pos=1)
}


VL_violin_boxplot <- function(w=100,h=100,loop_df,filepath, pair = "all"){
  if(pair != "all"){loop_df <- loop_df %>% filter(NAME == pair)}
  plots=list()
  for (i in 1:nrow(loop_df)) {
    HLA <- loop_df$hla[i]
    VAR <- loop_df$var[i]
    GENE <- loop_df$gene[i]
    NAME <- str_remove(loop_df$var[i], "x")
    colbox <- loop_df$col[i]
    hlaSNPbin <- fread(paste0("./hlaSNPbin_named_", GENE, ".csv")) %>% # generated in 03_hlaSNPbin_AA.R (called by 11_data_load.R)
      dplyr::select(ID, hla=all_of(HLA),var=all_of(VAR)) %>% drop_na(var) %>% dplyr::mutate(hla_var = as.factor(paste0(hla,"/",var))) # 00, 01, 10, 11

    hlaSNPvl <- full_join(hlaSNPbin,phenotype_B_naive%>%dplyr::select(ID,mean_logRNA,SEX), by="ID") %>% tibble() %>% drop_na(mean_logRNA,hla_var)
    
    give.n <- function(x){
      return(c(y = 7.5, label = length(x))) # experiment with the multiplier to find the perfect position
    }
    
    plot <- hlaSNPvl %>% drop_na() %>%
      ggplot(aes(y=mean_logRNA, x=hla_var))+
      geom_violin(width=0.8, fill="grey90", alpha=0.3,na.rm=T, color="grey40") +
      theme_classic() +
      
      labs(y=expression(paste(bold(~Log ["10"]~"viral load"))), 
           x= paste0(insert_colon(str_replace(HLA, "_", "*"))," / ",GENE, NAME))+
      stat_summary(fun.data = give.n, geom = "text", fun = median, size=4) +
      scale_y_continuous(breaks=c(0.0,2.0,4.0,6.0),labels=c("0.0","2.0","4.0","6.0"))+
      geom_jitter(position=position_jitter(0.4), size=0.5, col="grey50", alpha=0.6)+
      theme(axis.ticks.length = unit(2, units = "mm"),
            axis.text = element_text(size=12),
            panel.grid.major.y = element_line(colour = "grey90"),
            axis.title =  element_text(size=13, face = "bold"))+
      geom_boxplot(width=0.2, outlier.size = 0.5, col=colbox)
    
    cat(GENE, str_replace(HLA, "_", "*"),NAME,"\n")
    plots[[i]] <- plot
  }
  return(plots)
}

# Function to insert ":" at the second to last position
insert_colon <- function(x) {
  stringi::stri_sub(x, nchar(x)-1, nchar(x)-2) <- ":"
  return(x)
}

'%!in%' <- purrr::negate(`%in%`)

firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

roundUp <- function(x,m){
  if(x>0){return(m*ceiling(x / m))}
  if(x<0){return(m*floor(x / m))}
}

## alignement
source("./01_sequences_for_alignment.R")



