#!/usr/bin/env Rscript

# extract arguments
args <- commandArgs(trailingOnly = TRUE)

# Define arguments for testing
# [1]: where we get the isotype reference strains from sample_sheet.txt or CaeNDR release code, e.g. 20220216
# [2]: directory path for .bam and index files, default is "/projects/b1059/data/c_elegans/WI/alignments"
# [3]: path for the reference, default is "/projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
# [4]: path to R 3.6.0 libary
# args <- c("/projects/b1059/projects/Tim/sv-nf/temp_files/sample_sheet.txt",
#          "/projects/b1059/data/c_elegans/WI/alignments",
#          "/projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa")

# set libPath to AndersenLab R3.6.0 shared library for now
.libPaths(c(args[4], .libPaths() ))

library(magrittr)
#==============================================================================#
# Check how we should find reference isotypes                                 #
#==============================================================================#
if(!file.exists(args[1]) & !stringr::str_detect(args[1], pattern = "[0-9]{8}$")){
  
  # give a message and stop
  stop(glue::glue("sample_sheet.txt NOT found at {args[1]} - make sure the file is there and named sample_sheet.txt
                     or use --release argument with proper 8-digit release code"))
}

# detect sample sheet
if(stringr::str_detect(args[1], pattern = "sample_sheet.txt") & file.exists(args[1])){
  
  # give a message
  message(glue::glue("sample_sheet.txt found at {args[1]}"))
  
  # sort samples
  samples <- data.table::fread(args[1], header = F) %>%
    dplyr::rename(isotype = V1) %>%
    dplyr::arrange(isotype)
}

# detect release code
if(stringr::str_detect(args[1], pattern = "[0-9]{8}$")){
  
  # give a message
  message(glue::glue("pulling isotype reference strains from WI info sheet {args[1]}"))
  
  # get the species sheet
  samples <- gsheet::gsheet2tbl(url = "https://docs.google.com/spreadsheets/d/10x-CcKNCl80F9hMcrGWC4fhP_cbekSzi5_IYBY2UqCc") %>%
    dplyr::filter(release <= as.numeric(args[1])) %>%
    dplyr::filter(strain == isotype) %>%
    dplyr::arrange(isotype) %>%
    dplyr::select(isotype)
}

#==============================================================================#
# make the config file for delly                                               #
#==============================================================================#
config <- samples %>%
  dplyr::mutate(bam = paste0(args[2], "/", isotype, ".bam"),
                index = paste0(args[2], "/",  isotype, ".bam.bai"),
                ref = paste0(args[3]))

# write the config file
write.table(config, file = "delly_config_file.tsv", quote=FALSE, sep='\t', row.names = F)

# write the sorted sample file
write.table(samples, file = "sample_file.txt", quote=FALSE, sep='\t', row.names = F)

# reset the libPath - not sure if necessary
.libPaths(.libPaths()[-1])