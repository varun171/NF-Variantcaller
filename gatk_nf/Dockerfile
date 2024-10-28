# Base image - use an official Ubuntu image
FROM ubuntu:20.04

# Set the working directory inside the container
WORKDIR /pipeline

# Environment settings: Ensure non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    unzip \
    default-jre \
    build-essential \
    python3-pip \
    python3-dev \
    python-is-python3 \
    vim \
    samtools \
    bcftools \
    bwa \
    htop \
    && apt-get clean

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# Install GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && \
    unzip gatk-4.4.0.0.zip && \
    mv gatk-4.4.0.0 /usr/local/gatk && \
    ln -s /usr/local/gatk/gatk /usr/local/bin/gatk

# Install Python libraries for data manipulation and plotting (if needed)
RUN pip3 install --upgrade pip && \
    pip3 install pandas numpy matplotlib

# Copy Nextflow script and config file to the Docker container
COPY main.nf /pipeline/
COPY nextflow.config /pipeline/

# Add a simple healthcheck to verify that the key tools are installed properly
HEALTHCHECK --interval=5m --timeout=3s \
  CMD [ "/bin/bash", "-c", "gatk --help && bwa > /dev/null && nextflow -v && samtools --version" ] || exit 1

# Add entrypoint for easy access to Nextflow
ENTRYPOINT ["nextflow"]

# Expose port for Nextflow if you plan to use it with a web server
EXPOSE 8080

# Set the default command for the container
CMD ["-help"]
