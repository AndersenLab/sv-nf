FROM continuumio/miniconda
MAINTAINER Tim Crombie <tim.crombie@gmail.com>

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('IRanges')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('Rsamtools')"