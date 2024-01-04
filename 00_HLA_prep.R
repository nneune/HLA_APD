rm(list = ls())
source("packages_load.R") #packages

# HLA data ####
setwd("")
# call HLA data
...

##generates HLA 2 digit alleles from the 4 digit alleles
# for (i in 1:(ncol(hladata %>% dplyr::select(-COHORT,-freq,-rs333_T,-starts_with("PC"))))) {
#   # hladata[,i] <- str_pad(hladata[,i], 4, pad="0")
#   hladata[,i] <- substr(hladata[,i], 1,4)
# }

# remove letters
hladata <- hladata %>% ungroup() %>%
  dplyr::mutate(across(2:17, ~ str_replace_all(.,"[N|Q]","")))

# binomial
hlaBin <- hladata %>%
  dplyr::select(-COHORT, -rs333_T, -starts_with("PC")) %>%
  pivot_longer(cols = -ID, names_to = "variable", values_to = "value") %>%
  drop_na() %>% 
  dplyr::mutate(variable = ifelse(str_detect(variable, "_"), 
                                  str_remove(variable, "_.*"), 
                                  str_remove(variable, "\\d+")),
         cmb = paste(variable, value, sep = "_")) %>% 
  group_by(ID, cmb) %>%
  dplyr::summarise(dm = 1) %>% # ignore homolozygous alleles (=2)
  dplyr::arrange(cmb) %>%
  pivot_wider(names_from = cmb, values_from = dm, values_fill = 0L) %>%
  # mutate_all(~ ifelse(. > 0, 1, 0)) %>%
  replace(is.na(.), 0L) %>% 
  dplyr::arrange(ID)

# adding PCs again
hlaBin <- full_join(hlaBin, hladata %>% dplyr::select(ID,starts_with("PC"),COHORT,rs333_T), by="ID") 
fwrite(hlaBin, "hlaBin_4dig_corrected.csv",na = "NA")


# number of alleles ####
for (i in seq(2,16,by=2)) {
  nr_alleles <- unique(c(hladata[i]),c(hladata[i+1])) %>% table() %>% length()
  cat(colnames(hladata[i]), "has", nr_alleles, "alleles\n")
}
