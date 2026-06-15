FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="PARAFISH Platform"
LABEL description="Metabarcoding pipeline for fish parasites"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    unzip \
    build-essential \
    default-jre-headless \
    fasttree \
    libarchive-tools \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml
RUN conda env update --name base --file /tmp/environment.yml \
    && conda clean -afy

ARG VSEARCH_VERSION=2.28.1
RUN wget -q https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VERSION}/vsearch-${VSEARCH_VERSION}-linux-x86_64-static.tar.gz \
    && tar xzf vsearch-${VSEARCH_VERSION}-linux-x86_64-static.tar.gz \
    && mv vsearch-${VSEARCH_VERSION}-linux-x86_64-static/bin/vsearch /usr/local/bin/vsearch \
    && chmod +x /usr/local/bin/vsearch \
    && rm -rf vsearch-${VSEARCH_VERSION}*

# PEAR is installed via conda (see environment.yml). No binary download here.
# FastTree is installed via apt above; ensure the binary is available on PATH
RUN ln -sf $(which FastTree 2>/dev/null || which fasttree 2>/dev/null || echo /usr/bin/fasttree) /usr/local/bin/fasttree || true

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /app
COPY . /app

RUN mkdir -p /data/silva /data/pr2 /data/bold_db /data/blast_db /data/runs /data/config

# If conda installed Trimmomatic, create a stable path expected by the app
RUN mkdir -p /opt/Trimmomatic \
        && if [ -f /opt/conda/share/trimmomatic/trimmomatic.jar ]; then \
                 ln -sf /opt/conda/share/trimmomatic/trimmomatic.jar /opt/Trimmomatic/trimmomatic.jar; \
             fi \
        && if [ -d /opt/conda/share/trimmomatic/adapters ]; then \
                 ln -sf /opt/conda/share/trimmomatic/adapters /opt/Trimmomatic/adapters; \
             fi

ENV PATH="/opt/conda/bin:/usr/local/bin:$PATH"

EXPOSE 8050
CMD ["python", "app.py"]
