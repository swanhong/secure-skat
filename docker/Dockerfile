FROM gcr.io/google.com/cloudsdktool/cloud-sdk:slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      python3-numpy python3-scipy python3-matplotlib \
 && rm -rf /var/lib/apt/lists/*

COPY plink2 /usr/local/bin/plink2
RUN chmod +x /usr/local/bin/plink2 \
 && plink2 --version \
 && python3 -c 'import numpy, scipy, matplotlib'
