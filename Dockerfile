FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="PARAFISH Platform"
LABEL description="Metabarcoding pipeline for fish parasites"

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
# Dash listens on 127.0.0.1 unless HOST is set, which is unreachable from outside the container.
ENV HOST=0.0.0.0
ENV DEBIAN_FRONTEND=noninteractive
# The image has no .git folder; the commit is passed at build time for the run record
# (docker build --build-arg PARAFISH_COMMIT=$(git rev-parse HEAD) .).
ARG PARAFISH_COMMIT=unknown
ENV PARAFISH_COMMIT=$PARAFISH_COMMIT

# The base image is Debian 11 (bullseye). Its LTS ended in August 2026 and its packages
# moved to archive.debian.org; the regular mirrors now answer 404. bullseye-security was
# not archived, so it is dropped.
RUN printf 'deb http://archive.debian.org/debian bullseye main\ndeb http://archive.debian.org/debian bullseye-updates main\n' > /etc/apt/sources.list \
    && rm -f /etc/apt/sources.list.d/*.list /etc/apt/sources.list.d/*.sources \
    && apt-get -o Acquire::Check-Valid-Until=false update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    unzip \
    build-essential \
    default-jre-headless \
    fasttree \
    libarchive-tools \
    && rm -rf /var/lib/apt/lists/*

# The tools live in their own environment rather than in base: updating base would also
# have to update conda itself, which the solver may refuse. Strict channel priority keeps
# conda-forge and bioconda from mixing builds of the same library.
COPY environment.yml /tmp/environment.yml
RUN conda config --system --set channel_priority strict \
    && conda env create --name parafish --file /tmp/environment.yml \
    && conda clean -afy
ENV PATH="/opt/conda/envs/parafish/bin:/opt/conda/bin:/usr/local/bin:$PATH"

# VSEARCH comes from the Conda environment (environment.yml).

# PEAR is installed via conda (see environment.yml). No binary download here.
# FastTree is installed via apt above; ensure the binary is available on PATH
RUN ln -sf $(which FastTree 2>/dev/null || which fasttree 2>/dev/null || echo /usr/bin/fasttree) /usr/local/bin/fasttree || true

COPY requirements.txt /tmp/requirements.txt
# lxml_html_clean: lxml >= 6 moved lxml.html.clean out, and BOLDigger3 still imports it.
RUN pip install --no-cache-dir -r /tmp/requirements.txt lxml_html_clean

WORKDIR /app
COPY . /app

RUN mkdir -p /data/silva /data/pr2 /data/bold_db /data/blast_db /data/runs /data/config

# Stable paths for the Trimmomatic JAR and adapters expected by core/config.py. Bioconda
# installs them under a versioned folder, so it is located with a glob.
RUN mkdir -p /opt/Trimmomatic \
    && J=$(ls /opt/conda/envs/parafish/share/trimmomatic*/trimmomatic*.jar 2>/dev/null | head -1) \
    && if [ -n "$J" ]; then \
           ln -sf "$J" /opt/Trimmomatic/trimmomatic.jar; \
           ln -sfn "$(dirname "$J")/adapters" /opt/Trimmomatic/adapters; \
       else echo "WARNING: Trimmomatic JAR not found in the Conda environment"; fi


EXPOSE 8050
CMD ["python", "app.py"]
