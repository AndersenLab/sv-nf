#!/usr/bin/env Rscript

# extract arguments
args <- commandArgs(trailingOnly = TRUE)

# Define arguments for testing
# [1]: species, one of elegans, briggsae, tropicalis
# [2]: where we get the isotype reference strains from sample_sheet.txt or CaeNDR release code, e.g. 20220216
# [3]: directory path for .bam and index files, default is "/projects/b1059/data/c_elegans/WI/alignments"
# [4]: path for the reference, default is "/projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
# args <- c("elegans", "/projects/b1059/projects/Tim/sv-nf/temp_files/sample_sheet.txt",
#          "/projects/b1059/data/c_elegans/WI/alignments",
#          "/projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa")

library(magrittr)
#==============================================================================#
# Check how we should find reference isotypes                                 #
#==============================================================================#
if(!file.exists(args[2]) & !stringr::str_detect(args[2], pattern = "[0-9]{8}$")){
  
  # give a message and stop
  stop(glue::glue("sample_sheet.txt NOT found at {args[2]} - make sure the file is there and named sample_sheet.txt
                     or use --release argument with proper 8-digit release code"))
}

# detect sample sheet
if(stringr::str_detect(args[2], pattern = "sample_sheet.txt") & file.exists(args[2])){
  
  # give a message
  message(glue::glue("sample_sheet.txt found at {args[2]}"))
  
  # sort samples
  samples <- data.table::fread(args[2], header = F) %>%
    dplyr::rename(strain = V1) %>%
    dplyr::arrange(strain)
}

# detect release code
if(stringr::str_detect(args[2], pattern = "[0-9]{8}$")){
  # choose the species sheet to use
  if(args[1] == "elegans"){
    s <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/10x-CcKNCl80F9hMcrGWC4fhP_cbekSzi5_IYBY2UqCc")
  }
  if(args[1] == "briggsae"){
    s <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/1IJHMLwuaxS_sEO31TyK5NLxPX7_qSd0bHNKverAv8-0")
  }
  if(args[1] == "tropicalis"){
    s <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4")
  }
  if(args[1] != "elegans" & args[1] != "briggsae" & args[1] != "tropicalis"){
    # give a message
    stop(glue::glue("--species argument not recognized. Please set to one of elegans, briggsae, or tropicalis"))  
  }
  # give a message
  message(glue::glue("pulling isotype reference strains from {args[1]} WI info sheet"))
  
  # get the species sheet
  samples <- s %>%
    dplyr::filter(release <= as.numeric(args[2]) & isotype_ref_strain == T) %>%
    dplyr::arrange(strain) %>%
    dplyr::select(strain)
}

#==============================================================================#
# make the config file for delly                                               #
#==============================================================================#
config <- samples %>%
  dplyr::mutate(bam = paste0(args[3], "/", strain, ".bam"),
                index = paste0(args[3], "/",  strain, ".bam.bai"),
                ref = paste0(args[4]))

# write the config file
write.table(config, file = "delly_config_file.tsv", quote=FALSE, sep='\t', row.names = F)

# write the sorted sample file
write.table(samples, file = "sample_file.txt", quote=FALSE, sep='\t', row.names = F)