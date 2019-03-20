FROM continuumio/miniconda:4.5.4
LABEL authors="phil palmer phil@lifebit.ai" \
      description="Docker image containing the requirements for the svim pipeline"

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps gcc && apt-get clean -y 
RUN conda install conda=4.6.7

RUN conda config --add channels bioconda && \
    conda install -c bioconda svim