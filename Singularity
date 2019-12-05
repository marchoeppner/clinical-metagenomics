Bootstrap:docker
From:nfcore/base

%labels
    MAINTAINER Marc Hoeppner <m.hoeppner@ikmb.uni-kiel.de>
    DESCRIPTION Singularity image containing all requirements for the calinical metagenomics pipeline
    VERSION 1.0

%environment
    PATH=/opt/conda/envs/clinical-metagenomics-1.0/bin:$PATH
    export PATH

%files
    environment.yml /

%post

    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a

    apt-get -y update
    apt-get -y install procps build-essential

    apt-get -y install ruby ruby-dev imagemagick libmagickwand-dev
    gem install thinreports gruff

    mkdir -p /ifs

    mkdir -p /opt/trimmomatic && cd /opt/trimmomatic && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip \
        && unzip Trimmomatic-0.36.zip && mv Trimmomatic-0.36 0.36 && rm Trimmomatic-0.36.zip
