#!/usr/bin/env Rscript

# extract arguments
args <- commandArgs(trailingOnly = TRUE)

# Define arguments for testing
# [1]: bed file from delly /projects/b1059/projects/Tim/sv-nf/temp_files/WI.DELLYpif.germline.bed
# [2]: ref: /projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa
# [3]: path to R 3.6.0 libary
# args <- c("/projects/b1059/projects/Tim/sv-nf/temp_files/WI.DELLYpif.germline.bed", "/projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa", "/projects/b1059/software/R_lib_3.6.0")

# set libPath to AndersenLab R3.6.0 shared library for now
.libPaths(c(args[3], .libPaths() ))

# load packages
library(data.table)
library(Rsamtools)
library(magrittr)


# read in DELLY germline bed file
df <- data.table::fread(args[1]) %>%
  dplyr::select(CHROM = V1,
                START = V2,
                END = V3,
                SVTYPE = V4,
                STRAIN = V5,
                GT = V6) %>%
  dplyr::mutate(SIZE = END - START) %>%
  dplyr::filter(GT == "1/1")

# SAVE FORMATTED BED FILE
data.table::fwrite(df, "caendr.pif.bed", col.names = FALSE, sep="\t")
err = system("bgzip -f caendr.pif.bed && tabix -f caendr.pif.bed.gz")
if (err != 0) {
  Error("Something went wrong!")
}

STRAIN_ORDER <- sort(unique(df$STRAIN))

# Subset columns to only those that are needed.
vcf <- df[, c("CHROM", "START", "END", "SVTYPE", "SIZE", "STRAIN", "GT")]
vcf <- unique(vcf)
vcf <- maditr::dcast(vcf, CHROM + START + END + SVTYPE + SIZE ~ STRAIN, value.var="GT") # TAC maditr, could use dplyr

# Get reference sequences for ranges
REFERENCE <- glue::glue("{args[2]}") 
ref <- Rsamtools::scanFa(REFERENCE)
idx <- Rsamtools::scanFaIndex(REFERENCE)
seq_lengths <- IRanges::width(IRanges::ranges(idx))
names(seq_lengths) <- seqnames(idx)
range_set <- IRanges::IRanges(start = vcf$START, end = vcf$END)
ranges <- GenomicRanges::GRanges(seqnames = as.character(vcf$CHROM), ranges = range_set, seqlengths = seq_lengths)
ranges <- GenomicRanges::trim(ranges)
seqs <- Rsamtools::scanFa(REFERENCE, ranges)

# Integrate sequences
vcf[, sequence := as.character(seqs)]

# Setup VCF columns
vcf[, POS := START]
vcf[, ID := "."]

vcf[, REF := dplyr::case_when(SVTYPE == "INS" ~ substr(sequence, 1, 1),
                              SVTYPE == "DEL" ~ sequence), by=1:nrow(vcf)]
# insert size is - 1 of size for ref.
vcf[, ALT := dplyr::case_when(SVTYPE == "DEL" ~ substr(sequence, 1, 1),
                              SVTYPE == "INS" ~ paste0(substr(sequence, 1, 1), paste0(rep("A", SIZE), collapse=""), collapse="")), by=1:nrow(vcf)]
vcf[, QUAL := 1]
vcf[, INFO := paste0("INDEL=1;", "TYPE=", SVTYPE)]
vcf[, FILTER := "PASS"]
vcf[, FORMAT := "GT"]

vcf[, SVTYPE := NULL]
vcf[, SIZE := NULL]
vcf[, START := NULL]
vcf[, END := NULL]
vcf[, sequence := NULL]



setcolorder(vcf, c("CHROM",
                   "POS",
                   "ID",
                   "REF",
                   "ALT",
                   "QUAL",
                   "FILTER",
                   "INFO",
                   "FORMAT",
                   STRAIN_ORDER))
setnames(vcf, "CHROM", "#CHROM")


# Fix Genotypes; Set NA to 0/0; Infer reference
vcf[, (STRAIN_ORDER) := lapply(.SD, function(x) ifelse(is.na(x), "0/0", x)), .SDcols=STRAIN_ORDER]


HEADER_LINES <- c("##fileformat=VCFv4.2",
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##INFO=<ID=INDEL,Number=1,Type=Flag,Description="1 if indel">',
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="end position">',
                  '##INFO=<ID=TYPE,Number=1,Type=STRING,Description="type of variant">')

writeLines(HEADER_LINES, file("caendr.pif.vcf"))
data.table::fwrite(vcf, "caendr.pif.vcf", col.names=TRUE, append=TRUE, sep="\t")
err = system("bgzip -f caendr.pif.vcf && bcftools index caendr.pif.vcf.gz")
if (err != 0) {
  Error("Something went wrong!")
}

# reset the libPath - not sure if necessary
.libPaths(.libPaths()[-1])
