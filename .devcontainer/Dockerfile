FROM mcr.microsoft.com/vscode/devcontainers/base

# Install dependencies
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    # Common useful utilities
    gfortran \
    libncurses-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install CDF
ARG CDF_VER=39_0
RUN wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest/linux/cdf${CDF_VER}-dist-cdf.tar.gz && \
    tar -xvf cdf${CDF_VER}-dist-cdf.tar.gz && \
    cd cdf${CDF_VER}-dist && \
    make OS=linux ENV=gnu all -j && \
    make OS=linux ENV=gnu INSTALLDIR=/usr/local/cdf install && \
    cd .. && \
    rm -rf cdf${CDF_VER}-dist && \
    rm cdf${CDF_VER}-dist-cdf.tar.gz

ENV CDF_LIB /usr/local/cdf/lib
ENV SPEDAS_DATA_DIR /data