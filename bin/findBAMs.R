#!/usr/bin/env Rscript
library(tidyverse)

# extract arguments
args <- commandArgs(trailingOnly = TRUE)

# Define arguments for testing
# [1]: species, e.g. elegans, briggsae, tropicalis
# [2]: eight digit CenDR release code, e.g. 20220216
# [3]: name of output, e.g. 20220923_delly_indels.bcf
# [4]: type of sv to call, e.g. INS, DEL
#args <- c("elegans", "20220216", "20221108_delly_indels.bcf", "INS")

# pull species sheet
if(!(args[1] %in% c("elegans", "briggsae", "tropicalis"))){
  message("please choose a valid species: 'briggsae', 'elegans', or 'tropicalis'")
}
if(args[1] == "elegans"){
ss <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/10x-CcKNCl80F9hMcrGWC4fhP_cbekSzi5_IYBY2UqCc")
}
if(args[1] == "birggase"){
  ss <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/1IJHMLwuaxS_sEO31TyK5NLxPX7_qSd0bHNKverAv8-0")  
}
if(args[1] == "tropicalis"){
  ss <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4")  
}

# get isotype reference strains based on release
isos <- ss %>%
  dplyr::filter(release <= args[2]) %>%
  dplyr::filter(strain == isotype) %>%
  dplyr::mutate(isotype = paste0("/projects/b1059/data/c_elegans/WI/alignments/", isotype, ".bam")) %>%
  dplyr::pull(isotype)

# make a command to run delly, not sure what to do about home directory right now. This will change with user. NExtflow is the answer ultimately.
# this one will run independent commands - could be useful to plug into a nextflow channel.
#glue::glue("singularity exec -B /projects/b1059/projects/Tim/sv-nf:/home/tac9384 delly_latest.sif delly call -t {args[4]} -o {args[3]} -g temp_files/c_elegans.PRJNA13758.WS283.genome.fa {isos}")

# make the command
com_INS <- glue::glue("singularity exec -B /projects/b1059/data/c_elegans/WI/alignments,/projects/b1059/projects/Tim/sv-nf:/home/tac9384 delly_latest.sif delly_latest.sif delly call -t {args[4]} -o {args[3]} -g temp_files/c_elegans.PRJNA13758.WS283.genome.fa {paste(isos, collapse = ' ')}")
com_DEL <- glue::glue("singularity exec -B /projects/b1059/data/c_elegans/WI/alignments,/projects/b1059/projects/Tim/sv-nf:/home/tac9384 delly_latest.sif delly call -t DEL -o {args[3]} -g temp_files/c_elegans.PRJNA13758.WS283.genome.fa {paste(isos, collapse = ' ')}")

# save the command to run
write_file(com_INS, file = "/projects/b1059/projects/Tim/nil-ril-nf/bin/com_INS.txt")
write_file(com_DEL, file = "/projects/b1059/projects/Tim/nil-ril-nf/bin/com_DEL.txt")
