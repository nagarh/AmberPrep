FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    build-essential \
    gcc \
    g++ \
    make \
    libffi-dev \
    libssl-dev \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:${PATH}"

# Accept conda Terms of Service and configure channels
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority flexible

# Install mamba for faster package installation
RUN conda install -n base -c conda-forge mamba -y

# Install AMBER tools and PyMOL using mamba with Python 3.11 (compatible version)
RUN mamba install -y python=3.11 conda-forge::ambertools conda-forge::pymol-open-source

# Clean up conda/mamba cache to reduce image size
RUN conda clean -afy

# Install Python packages via pip
RUN pip install --no-cache-dir \
    flask==2.3.3 \
    flask-cors==4.0.0 \
    biopython==1.81 \
    numpy==1.24.3 \
    pandas==2.0.3 \
    matplotlib==3.7.2 \
    seaborn==0.12.2 \
    mdanalysis==2.5.0 \
    gunicorn==21.2.0 \
    requests==2.31.0 \
    rdkit==2023.3.1 \
    scipy==1.11.1

# Set working directory
WORKDIR /AmberFlow

# Copy the entire project
COPY . .

# Create necessary directories with proper permissions
RUN mkdir -p /AmberFlow/obsolete /AmberFlow/pdb /AmberFlow/temp /AmberFlow/output && \
    chmod -R 777 /AmberFlow

# Make sure the python directory is in the Python path
ENV PYTHONPATH="${PYTHONPATH}:/AmberFlow/python"

# Expose the port
EXPOSE 7860

# Run the application
CMD ["python", "start_web_server.py"]

