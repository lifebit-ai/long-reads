FROM continuumio/miniconda:4.5.4
LABEL authors="phil palmer phil@lifebit.ai" \
      description="Docker image containing base requirements for the clairvoyante pipeline"

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps gcc && apt-get clean -y 
RUN conda install conda=4.6.7

RUN pip install tensorflow==1.9.0 && \
    pip install blosc && \
    pip install intervaltree==2.1.0 && \
    pip install numpy

RUN conda config --add channels conda-forge && \
    conda install -c conda-forge pypy2.7==5.10.0 && \
    conda install -c conda-forge python-blosc==1.8.1 && \
    conda install -c conda-forge intervaltree==2.1.0

RUN wget https://bootstrap.pypa.io/get-pip.py && \
    pypy get-pip.py && \
    pypy -m pip install --no-cache-dir intervaltree==2.1.0

RUN conda config --add channels bioconda && \
    conda install -c bioconda clairvoyante && \
    clairvoyante.py

RUN apt-get install parallel -y

RUN conda install -c bioconda samtools openssl=1.0 && \
    conda install -c bioconda htslib && \
    conda install -c bioconda vcflib