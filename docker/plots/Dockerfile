FROM continuumio/miniconda:4.5.4
LABEL authors="phil palmer phil@lifebit.ai" \
      description="Docker image containing the requirements for structural varaints plots"

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps gcc && apt-get clean -y 
RUN conda install conda=4.6.7

RUN conda config --add channels bioconda && \
    conda install -c bioconda cyvcf2

RUN conda config --add channels conda-forge && \
    conda install -c conda-forge matplotlib

COPY ./SV-carriers-plot.py /usr/local/bin
COPY ./SV-length-plot.py /usr/local/bin

RUN chmod +x /usr/local/bin/SV-carriers-plot.py && \
	chmod +x /usr/local/bin/SV-length-plot.py

ENV PATH=$PATH:/usr/local/bin

