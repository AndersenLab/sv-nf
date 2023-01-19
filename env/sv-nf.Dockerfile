FROM continuumio/miniconda
MAINTAINER Tim Crombie <tim.crombie@gmail.com>

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

# add command line utilities to image
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps