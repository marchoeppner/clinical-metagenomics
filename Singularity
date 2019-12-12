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

    apt-get -y install ruby ruby-dev imagemagick libmagickwand-dev libncurses5-dev
    gem install thinreports gruff

    mkdir -p /ifs

    cd /opt
    mkdir -p samtools
    cd samtools
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar -xvf samtools-1.9.tar.bz2
    mv samtools-1.9 source 
    cd source
    ./configure --prefix=/opt/samtools/1.9 && make install
    cd /opt/samtools
    rm -Rf source *.bz2
