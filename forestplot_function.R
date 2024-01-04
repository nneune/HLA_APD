forestploter_table <- function( # structure of data: HLA,VAR,gene,term,estimate,interaction_sign,conf.low,conf.high,pval,nobs,orCI,POS,hla,var,var_names
  data_name = "df_name", 
  ref_all = T, 
  gene, 
  plotname=NULL,
  legend.pos = list("top"),
  folder = paste0("path"),
  size=1,
  batch_size = 8, # Adjust the batch size as needed,
  gap=8,
  lineheight = "auto",
  graphwidth = "auto",
  # return=F,
  line.margin=NULL,
  show_legend=T,
  g_width = 1650,
  g_height=1300,
  graph.pos=3
){
  if(!grepl(data_name, pattern = "VL|APD")){stop("wrong model input")}
  ma <- get(data_name) %>% dplyr::select("HLA","VAR","region"="gene","term","estimate","interaction_sign","nobs","conf.low","conf.high","pval","orCI","POS", "hla", "var", "var_names")
  if(gene != "all"){ma <- ma %>% filter(region==gene) %>% ungroup()}else{ma <- ma %>% ungroup()}
  if (!file.exists(folder)) {
    dir.create(folder)
  } 
  
  ma <- ma %>% mutate(term=as.factor(term))

  if(nrow(ma)==0){return(NULL)}
  ok <- TRUE
  k <- 0

  while (ok) {
    if (nrow(ma) / (batch_size*4) > 1) {
      start_row <- k * (batch_size*4) + 1
      end_row <- (k + 1) * (batch_size*4)
      
      m <- ma %>%
        filter(row_number() >= start_row & row_number() <= end_row) %>%
        group_by(HLA, VAR) %>%
        arrange(term, .by_group = TRUE)
      
      if (nrow(m) < 4) {
        break
      }
      cat(gene, " ")
      k <- k + 1
    } else {
      m <- ma %>%
        group_by(HLA, VAR) %>%
        arrange(term, .by_group = TRUE) 
      cat(gene, " ")
      ok <- FALSE
    }
    
    m <- m %>% arrange(region,POS)
    
    if(grepl(data_name, pattern = "APD")){
      expo = T
      main = "Amino Acid association with HLA allele - APD interaction"
      nrows_forest = 4
      legend_text <- as.matrix(m %>% ungroup() %>%dplyr::mutate(region= as.character.factor(region)) %>%
                                 dplyr::reframe(legend = paste0(ifelse(gene=="all", paste0(firstup(region),var), var),
                                                                " ~ ", hla,str_pad(width = (20-(nchar(VAR)+nchar(HLA))),
                                                                                   string = paste0("(n=",nobs,")"), side = "left")),.by = VAR) %>% dplyr::select(-VAR) %>% distinct())
      xt=log(c(1e-1,0.5,1,2,5,2e1,1e2))
      attr(xt, "labels") <- as.character(exp(xt))

    }
    if(grepl(data_name, pattern = "VL")){
      nrows_forest = 3
      expo = F
      main = "Viral load association with HLA allele - HIV variant interaction"
      legend_text <- as.matrix(unique(m %>% ungroup() %>% arrange(term) %>% dplyr::mutate(region= as.character.factor(region)) %>%
                                        dplyr::reframe(legend = paste0("VL ~ ", paste0(region,var)," x ",hla,str_pad(width = (20-(nchar(VAR)+nchar(HLA))),string = paste0("(n=",nobs,")"), side = "left")))))
      xt=seq(-1.5, 1.5, 0.5)
      attr(xt, "labels") <- as.character(xt)
    } 
    
    m <- m %>% 
      {if (grepl(data_name, pattern = "VL")) filter(.,term!="APD") else .}
    
    bindOR = m %>% 
      ungroup() %>% arrange(term) %>% 
      dplyr::mutate(title = ifelse(term == lag(term) & row_number()>1, "",as.character.factor(term))) %>%
      dplyr::mutate(var_names = ifelse(var_names == term, "",var_names)) %>% 
      
      dplyr::summarise(
        .by = term,
        cbindorCI = paste(orCI, collapse = "\n"),
        cbindpval = paste(pval, collapse = "\n"),
        names = paste(var_names, collapse = "\n"),
        title_names = paste(title, collapse = "\n")
      )
    
    if(grepl(data_name, pattern = "APD")){
      bindOR$title_names[2] <- paste0("APD if HLA is \n\n\n\n\n\n\n")
      bindOR$title_names[3] <- paste0("APD if HLA is \n\n\n\n\n\n\n")
      bindOR$names[2] <- paste0("absent\n\n\n\n\n\n\n")
      bindOR$names[3] <- paste0("present\n\n\n\n\n\n\n")}
    
    # Specifying plot description
    text= bindOR %>% dplyr::relocate(title_names,names,.before=term) %>% 
      dplyr::select(-cbindpval) %>%
      dplyr::select(-term)   %>%
      {if (show_legend==T) dplyr::mutate(.,`         `=str_c(rep(" ", nchar(legend_text[1])),collapse = "")) else .}
    
    estimates <- m %>% 
      ungroup() %>% arrange(term) %>% 
      summarise(
        .by = term,
        me= list(estimate),
        low= list(conf.low),
        up= list(conf.high)
      )
    
    if(show_legend==F){legend_text=NULL}
    gap <- as.numeric(gap)
    size <- as.numeric(size)
    graph.pos <- as.numeric(graph.pos)
    plot_save_as=plotname
    if(is.null(plotname)){plot_save_as = paste(gene,k,ifelse(expo==T, "OR", "beta"),ifelse(ref_all == T, "ref_all","ref_dot"),plotname,"forest.tiff", sep = "_")}
    
    tiff(filename = paste0(folder,plot_save_as), width = g_width, height = case_when(nrow(m)<(batch_size/2)*nrows_forest~g_height-400,nrow(m)>=(batch_size/2)*nrows_forest~g_height),units = "px", compression = "none",res = 100)
    x <- T
    tryCatch({
      fp <- 
        (forestplot::forestplot(
          text,
          labeltext=c(title_names,names,cbindorCI),
          mean = rbind(matrix(unlist(estimates$me),nrow = length(estimates$me), byrow = T)),
          lower = rbind(matrix(unlist(estimates$low),nrow = length(estimates$low), byrow = T)),
          upper = rbind(matrix(unlist(estimates$up),nrow = length(estimates$up), byrow = T)),
          align = c("l", "l","r","r"),# text right, left, center
          is.summary = F, # which text is bold which is not
          graph.pos = graph.pos,  # position of graph in between the data text
          clip=c(4e-2, 9e1),
          xticks.digits = 3, # number of decimals
          xlog = expo, # TRUE if OR
          xticks =xt,
          line.margin=line.margin,
          txt_gp = forestplot::fpTxtGp(cex = size, # text size
                                       label = gpar(fontfamily="Arial",cex=size),
                                       legend= grid::gpar(cex=size,line.margin=line.margin+0.5), # legend text size
                                       ticks = grid::gpar(cex=size), # number size at ticks
                                       xlab = grid::gpar(cex=size+0.2, fontface = "bold")), # xlab text size
          xlab = ifelse(expo==T,"Odds ratio of the HIV variant to occur","Effect on Viral Load"),
          
          legend = legend_text, # legend text
          legend_args = forestplot::fpLegend(
            title = "Variant pairs",
            pos = legend.pos,
            gp = grid::gpar(col = "#FFFFFF")),  #"#F9F9F9")),
          
          {if(show_legend==F)legend=NULL},
          lwd.xaxis = grid::gpar(lwd=1.6),
          lwd.zero = grid::gpar(lwd=1.6), # thickness of zero line
          lwd.ci = (size+0.4), 
          ci.vertices = TRUE, 
          ci.vertices.height = (size*0.02), # confidence interval size and thickness
          lineheight = lineheight,
          graphwidth = graphwidth,
          
          boxsize = (size*0.03), # size of data point
          colgap = unit(gap, "mm"), # distance between cols
          col = forestplot::fpColors(box = c("#332288","#0066CC", "#12A4EA","#117733","#EFB412","#C00000","#CC0066", "#731941"),
                                     lines = c("#332288","#0066CC", "#12A4EA","#117733","#EFB412","#C00000","#CC0066", "#731941"),
                                     zero = c("black", "#332288"))) %>% 
           fp_add_header(title_names= "\n\n\nVariable",
                         names = "\n ",
                         cbindorCI= ifelse(grepl(data_name, pattern = "APD"),"\n\n\nOR [95% CI]","\n\n\nβ [95% CI]") 
                         # ,cbindpval="P-value"
           ) 
        )
        
      }, error = function(e) {
          message(paste0(gene, "-"),conditionMessage(e))
          
          x <<- F
        })
    if(return){return(fp)}
    if(x){cat("Picture saved \n")}
    print(fp)
    invisible(capture.output(dev.off()))
  }
}
