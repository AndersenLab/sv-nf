library(tidyverse)
#BiocManager::install("intansv")
library(intransv)

# open the bed file
bed <- data.table::fread("/projects/b1059/projects/Tim/sv-nf/out/CeNDR_pairwise_indel_finder/sv.20200815.bed.gz") %>%
  dplyr::mutate(size = V3 - V2) %>%
  dplyr::group_by(V8) %>%
  dplyr::mutate(n.samples = n()) %>%
  dplyr::ungroup()

# get size range for INS and DEL
max(bed$V17)
max(bed$size)
min(bed$size)

# check unique SV types
types <- unique(bed$V5)
types <- table(bed$V7)

# see how many singletons there are: some, also many with mutliple records per sample because it looks 
ggplot(bed) +
  aes(x = n.samples) +
  geom_histogram()
#==============================================================================#
# Read in delly output file with intransv - NOT WORKING YET
#==============================================================================#
#del.out <- intansv::readDelly(file = "/projects/b1059/projects/Tim/sv-nf/out/test/germline.bcf", method = "Delly")
str(delly)
