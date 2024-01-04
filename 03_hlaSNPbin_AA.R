## called within 11_data_load.R ## region is replaced by respective HIV gene
cat("region")
binary_region_aa_alt <- fread(paste0(folder,"binary_region_aa_alt_dot")) ## generated in 02_binary_region_aa.R

SNP_region_table <- binary_region_aa_alt %>% remove_rownames() %>% group_by(POS) %>% filter(ALT!= "Z") %>% dplyr::select(-REF, -ALT) # remove all gaps in all Alternatives
SNP_region_table$rownumber = 1:nrow(SNP_region_table)
SNP_region_table <- SNP_region_table %>% dplyr::relocate(rownumber,.before = POS) %>% dplyr::rename(id = ID)
SNP_region_table <- SNP_region_table %>% dplyr::mutate(ID = paste0("V",rownumber)) %>%
  dplyr::relocate(ID,.before = POS) %>% group_by(ID) %>% dplyr::select(-rownumber)

# backup table as reference 
SNP_region_table_ref <- SNP_region_table
write.csv(SNP_region_table_ref, paste0(folder,"SNP_region_table_ref_dot.csv"))
cat("SNP_region_table_ref, your translation for V1-Vn, is saved to", getwd(),"\n")
SNP_region_table <- SNP_region_table %>% dplyr::select(-POS, -id)
SNP_region_table <- (t(SNP_region_table) %>% as.data.frame() %>% rownames_to_column("base_uuid"))[-1,]
write.table(SNP_region_table,paste0(folder,"SNP_region_table_dot.csv"), sep ="\t", quote = FALSE, row.names = FALSE)

# merging HLA and SNP table
if(!"base_uuid" %in% colnames(hlaBin)){hlaBin <- merge(hlaBin, NGS %>% select(ID, base_uuid))}
hlaSNPbin_region <- merge(SNP_region_table, hlaBin, by="base_uuid", all=F)
hlaSNPbin_region <- hlaSNPbin_region %>% relocate(ID, .before = base_uuid)
# for dot data
hlaSNPbin_region[hlaSNPbin_region == "."] = NA

hlaSNPbin_region[3:ncol(hlaSNPbin_region)] <- mutate_if(hlaSNPbin_region[3:ncol(hlaSNPbin_region)], is.character, as.numeric)
hlaSNPbin_region <- distinct(hlaSNPbin_region)

fwrite(na="NA",hlaSNPbin_region, "./hlaSNPbin_region_dot.csv")

cat("hlaSNPbin_region_dot.csv is saved to", getwd(),"\n")

if(keep_files == F){rm("hlaSNPbin_region","SNP_region_table","SNP_region_table_ref")}

cat("_______________________________________________________________________________________________________________________________________________________\n\n")

