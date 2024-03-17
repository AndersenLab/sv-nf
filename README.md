# sv-nf

*A nextflow pipeline to call structural variants in Caenorhabditis nematodes from .bam files*

## Pipeline overview
<p align="center" width="100%">
    <img width="60%" src="https://github.com/AndersenLab/sv-nf/blob/main/img/sv-nf_workflow.png?raw=true">
</p>

## QUEST usage

The current version of `sv-nf v0.1.0` is limited to building the `.bed` and `.vcf` files required for the CaeNDR pairwise-indel finder tool. Follow the usage below to generate the files. 

```
# clone the repo
git clone <https://github.com/AndersenLab/sv-nf.git>
cd sv-nf

# setup environment
module load python/anaconda3.6
module load singularity
source activate /projects/b1059/software/conda_envs/nf20_env/

# example run for latest CeNDR release
nextflow run main.nf --sepcies elegans --ref /projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa --release 20220216
```

## sv-nf help

```bash
S V - N F    P I P E L I N E
===============================================
Usage:
The typical usage is as follows:
nextflow run main.nf --sepcies elegans --ref /projects/b1059/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa --release <latest CaeNDR release>
To debug use:
nextflow run main.nf --debug

Required Arguments:
--species     String           One of three, elegans, briggsae, tropicalis
--ref         String           Full path to the .fa uncompressed reference file, default is null
--release     String           The 8-digit CaeNDR release date, e.g 20220216, required when --sp_sheet not specified
OR
--sp_sheet    String           A path to the sample_sheet.txt file used to call INDELs instead of release. Required when --release not specified

Optional Arguments:
--bam_dir     String           The path to the .bam directory, default is specific to QUEST: /projects/b1059/data/c_<species>/WI/alignments
--out         String           The output directory, default is SV_indel_results_<date>
--debug

Flags:
--help                          Display this message

Notes:
The briggsae reference path on QUEST is /projects/b1059/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome.fa.gz
The tropicalis reference path on QUEST is /projects/b1059/data/c_tropicalis/genomes/NIC58_nanopore/June2021/c_tropicalis.NIC58_nanopore.June2021.genome.fa.gz
```
