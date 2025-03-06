# Nick Sunderland March 2025
# nicholas.sunderland@bristol.ac.uk

# base image (ARM here as I am on Mac m2)
FROM arm64v8/python:3.7-bullseye
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    openjdk-11-jdk-headless \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm -rf /tmp/miniconda.sh && \
    /opt/conda/bin/conda init bash

# Set up JAVA environment (ARM here as I am on Mac m2)
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-arm64
ENV PATH="${JAVA_HOME}/bin:${PATH}"

# Set up Conda environment
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /usr/src/app

# Copy the environment.yml file to the container
COPY environment.yml .

# Create the Conda environment based on the environment.yml file
RUN /opt/conda/bin/conda env create -f environment.yml

# Copy the rest of the application files
COPY . .

# Define the entrypoint for the container
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "gwas_qc", "python", "main.py"]

# Build instructions
# To build:
# cd this_Dockerfile/dir/path
# docker build --progress=plain --tag nicksunderland/gwas_qc --file Dockerfile . &> docker_build.log

# push to dockerhub
# docker push nicksunderland/gwas_qc:latest