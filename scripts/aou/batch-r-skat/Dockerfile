FROM gcr.io/deeplearning-platform-release/base-cpu@sha256:394be02c6b020a39837e0719c546d56b829994ce091fde7970c232b7f16a6640

USER root
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        r-base r-base-dev build-essential gfortran \
    && Rscript -e 'install.packages("SKAT", repos="https://cloud.r-project.org")' \
    && Rscript -e 'library(SKAT)' \
    && rm -rf /var/lib/apt/lists/*
