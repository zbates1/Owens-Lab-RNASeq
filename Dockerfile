# Use rocker/tidyverse as base image
FROM rocker/r-base:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and required packages
RUN R -e 'install.packages("tidyverse", repos="https://cloud.r-project.org/")' && \
    R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org/")' && \
    R -e 'BiocManager::install(c("DESeq2", "pheatmap"), update = FALSE)'

# Create directories
RUN mkdir -p /01_data \
    /02_code \
    /03_output

# Copy files
COPY test.r /02_code/myScript.R
COPY /Col1GFP_Sepsis/genes.readcount.xls /01_data/genes.readcount.xls

# Set working directory
WORKDIR /02_code

# Set the entrypoint to run R script
ENTRYPOINT ["Rscript", "myScript.R"]
