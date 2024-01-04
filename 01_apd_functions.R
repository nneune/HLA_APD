##blast
region_extract = function(uuids,region) {
    region = tolower(region)
  ##BLAST database
  # rBLAST::makeblastdb(paste0("./references_seqs/hxb2_",region,".fa"),dbtype = "nucl") ## hxb2 ref
  # bl <- rBLAST::blast(paste0("./references_seqs/hxb2_",region,".fa"), type = "blastn")
  rBLAST::makeblastdb(paste0("./hxx_references/hxx_",region,".fa"),dbtype = "nucl") ## Los Alamos refs
  bl <- rBLAST::blast(paste0("./hxx_references/hxx_",region,".fa"), type = "blastn")
  blastmeta = data.frame(QueryID = character(),
                               SubjectID = character(),  
                               Perc.Ident = numeric(),  
                               Alignment.Length = integer(),
                               Mismatches = integer(),   
                               Gap.Openings = integer(),  
                               Q.start  = integer(),      
                               Q.end   = integer(),        
                               S.start  = integer(),        
                               S.end   = integer(),         
                               E = numeric(),             
                               Bits = numeric(), 
                               stringsAsFactors = FALSE)
  
  j=1
  for(i in uuids){
    cat("\n", i, "\n")
    ##sequences
    if(!file.exists(paste0(consensus_path,i,".fa"))){next} # check if file there otherwise skip
    # read seq (Biostrings package)
    seq <- Biostrings::readDNAStringSet(paste0(consensus_path,i,".fa"))
    seq_df = as.data.frame(Biostrings::readDNAStringSet(paste0(consensus_path,i,".fa")))
    seq <- Biostrings::chartr(seq, old = "-", new = "N")
    
    ##BLAST
    blast_seq = tryCatch({predict(bl, seq, BLAST_args = "-max_target_seqs 2000 -evalue 4000")},
                         error = function(e) {
                           print("why")
                         })
    if(nrow(blast_seq) == 0){
      blast_seq [ nrow(blast_seq) + 1 , ] <- NA
      blast_seq[1,1] = i
      blastmeta = rbind(blastmeta,blast_seq)
      
      cat(paste0(j, ' of ', length(uuids), ' completed'))
      j= j + 1
      next }
    
    blast_seq = blast_seq %>% dplyr::arrange(E,desc(Alignment.Length)) %>%
      dplyr::slice(1) # first with highest Alignment Length
    
    ##MERGE BLAST with raw sequence
    folder=paste0("./APD/blasted_seqs_",region,"/")
    if (!file.exists(folder)) {dir.create(folder)} # create folder for blasted seqs
    
    seq_df %>%
      mutate(reg_seq = substr(x, blast_seq$Q.start, blast_seq$Q.end))  %>%
      mutate(forcsv = paste0(">",rownames(seq_df)[1],"\n",reg_seq)) %>% {
      write.table(.$forcsv,paste0(folder,i,"_blasted.fa"), 
                  row.names = FALSE, quote = FALSE,col.names = FALSE)
      }
    blastmeta = rbind(blastmeta,blast_seq)
    cat(paste0(j, ' of ', length(uuids), ' completed\n'))
    j= j + 1
  }
  folder=paste0("./APD/blastmeta/")
  if (file.exists(folder)) {cat("\n")} else {dir.create(folder)}
  fwrite(blastmeta,paste0("./APD/blastmeta/blastmeta_",region,".csv"),na="NA")
}


##MACSE
condon_align = function(uuids,region){
  region = tolower(region)
  j=0
  folder=paste0(region,"_codon_align/")
  if (file.exists(folder)) {
    cat("\n")
  } else {dir.create(folder)}
  
  for(i in uuids){
    if (sum(file.exists("./APD/blasted_seqs_",region,"/",i,"_blasted.fa"))>0) {
      file = paste0("./APD/blasted_seqs_",region,"/",i,"_blasted.fa")
    } else {next}
    
    # MACSE run translation NT -> AA
    system(paste0("java -jar ",macse_path,
                  " -prog alignSequences -seq ./APD/references_seqs/hxb2_",region,".fa -seq_lr ", # ref is hxb2
                  file," -fs_lr 10 -stop_lr 15 -fs 100 -fs_term 100 ",
                  " -out_AA ./",folder,i,"_",region,"_AA.fa",
                  " -out_NT ./",folder,i,"_",region,"_NT.fa"))
    j=j+1
  }
}



##calculate score 
apd_calc = function(uuid,region,depth_threshold){
    region = tolower(region)

  ##count third codon position
  co_al = as.data.frame(read.table(paste0("./",region,"_codon_align/",uuid,"_",region,"_NT.fa"))[4,])
  colnames(co_al) = "seq"
  
  co_al = co_al %>% separate(seq, into = paste("POS", 0:(max(nchar(co_al$seq))), sep = ""), sep = "") %>% 
    dplyr::select(-'POS0') %>% tidyr::gather(number , base) %>% 
    dplyr::mutate(number =  rownames(.))
  
  co_al = co_al %>% 
    dplyr::filter(base %in% c("A","G","C","T","!", "N"))

  co_al = co_al %>% 
    dplyr::mutate(co_pos = rep(seq(1:3),nrow(.)/3)) %>%
    dplyr::mutate(pos_codon = paste0(number,"_",co_pos)) %>%
    dplyr::filter(base %in% c("A","G","C","T", "N"))

  ##get original assembly position
  seq_org_pos <- Biostrings::readDNAStringSet(paste0(consensus_path,uuid,".fa"))
  seq_meta = fread(paste0("./APD/blastmeta/blastmeta_",region,".csv")) %>% filter(grepl(QueryID, pattern=uuid))
  
  seq_org_pos = as.data.frame(chartr("N", "N", seq_org_pos))
  colnames(seq_org_pos) = "seq"
  seq_org_pos = seq_org_pos %>% 
    separate(seq, into = paste("POS", 0:(max(nchar(seq_org_pos$seq))), sep = ""), sep = "") %>% 
    dplyr::select(-'POS0') %>% tidyr::gather(number_assembly , base_assembly) %>% 
    dplyr::mutate(number_assembly =  rownames(.)) %>%
    dplyr::filter(number_assembly %in% seq_meta$Q.start:seq_meta$Q.end)
  
  ambg_pos = as.data.frame(cbind(co_al,seq_org_pos%>% filter(base_assembly != "-"))) %>% 
    filter(co_pos == 3) %>% filter(base != "N") # third position, no gaps

##load assembly information
  freqs = fread(paste0(freq_path,uuid,"_freq.csv"))
  covs = fread(paste0(cov_path,uuid,"_cov.csv"))
  covs = covs %>% filter(depth >= depth_threshold)

##calculate score
  ambg_score = freqs %>% filter(position_1basedindexing %in% ambg_pos$number_assembly) %>%
    dplyr::filter(position_1basedindexing %in% covs$position) %>%
    dplyr::group_by(position_1basedindexing) %>%
    dplyr::mutate(countit = max(frequency) < 0.99) %>%
    dplyr::mutate(score_base = frequency * (1-frequency)) %>%
    dplyr::mutate(ambg_score_pos = sum(score_base)) %>%
    dplyr::mutate(ambg_score_pos = ifelse(countit == FALSE,0,ambg_score_pos)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n = sum(countit),
                  apd_score = sum(ambg_score_pos)/n()) %>%
    dplyr::slice(1)

  return(ambg_score %>% select(apd_score,n)) # APD score and number of positions is was called over
}
