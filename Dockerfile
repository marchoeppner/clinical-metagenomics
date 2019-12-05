FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for clinical metagenomics pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/clinical-metagenomics-1.0/bin:$PATH

RUN apt-get -y update &&  apt-get -y install procps build-essential && apt-get -y install ruby ruby-dev imagemagick libmagickwand-dev && gem install thinreports gruff

RUN mkdir -p /opt/trimmomatic && cd /opt/trimmomatic && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip \
        && unzip Trimmomatic-0.36.zip && mv Trimmomatic-0.36 0.36 && rm Trimmomatic-0.36.zip
