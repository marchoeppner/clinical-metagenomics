FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for clinical metagenomics pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/clinical-metagenomics-1.0/bin:$PATH

RUN apt-get -y update &&  apt-get -y install procps build-essential && apt-get -y install ruby ruby-dev imagemagick libmagickwand-dev && gem install thinreports gruff

RUN cd /opt && \
    mkdir -p samtools && \
    cd samtools && \
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xvf samtools-1.9.tar.bz2 && \
    mv samtools-1.9 source && \
    cd source && \
    ./configure --prefix=/opt/samtools/1.9 && make install && \
    cd /opt/samtools && \
    rm -Rf source *.bz2
