FROM ubuntu:20.04
RUN apt update -y
RUN apt install -y software-properties-common
RUN add-apt-repository universe
RUN apt install -y python3.9 python3-pip
RUN apt install -y libjpeg-dev zlib1g-dev
RUN apt-get -y install libssl-dev libfontconfig1-dev libxml2-dev libcurl4-openssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libcairo2-dev pkg-config python3-dev openssl procps bcftools bedtools tabix r-base wget rsync
RUN pip3 install --upgrade pip setuptools wheel --default-timeout=100
RUN pip3 install pycairo
RUN pip3 install matplotlib
RUN pip3 install seaborn
RUN pip3 install SigProfilerMatrixGenerator
RUN pip3 install SigProfilerSimulator
RUN pip3 install SigProfilerClusters
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('Biostrings')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('grid')"
RUN R -e "install.packages('stringi')"
RUN R -e "install.packages('ggnewscale')"
RUN R -e "install.packages('purrr')"
RUN R -e "install.packages('stringr')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('VariantAnnotation')"
RUN R -e "BiocManager::install('StructuralVariantAnnotation')"
RUN R -e "BiocManager::install('annotate')"
RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "options(timeout = 6000); BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')"
RUN R -e "options(timeout = 6000); BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')"
RUN R -e "options(timeout = 6000); BiocManager::install('org.Hs.eg.db')"
RUN python3 -c "exec(\"from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh37', rsync=False)\")"
