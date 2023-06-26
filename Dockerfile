###############################################
# Dockerfile to build ClassifyCNV container image
#   docker build -t gvcn/classifycnv:v.1.1.1 . --platform linux/amd64
#   docker push classifycnv:v.1.1.1
#   docker run -v $PWD:$PWD classifycnv:v.1.1.1 python3 ClassifyCNV.py -h
#   mkdir output/; docker run -v $PWD:$PWD classifycnv:v.1.1.1 python3 ClassifyCNV.py --infile Examples/ACMG_examples.hg19.bed --GenomeBuild hg19 --precise --outdir $PWD/output
###############################################
FROM continuumio/miniconda3:23.3.1-0-alpine

# Install git
RUN conda install -c anaconda git && \
    conda install -c bioconda bedtools && \
    git clone --branch v1.1.1 --depth 1 https://github.com/Genotek/ClassifyCNV.git && \
    cd ClassifyCNV && \
    bash update_clingen.sh

WORKDIR "./ClassifyCNV"
ENV PATH ${PWD}:${PATH}

