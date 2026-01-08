FROM continuumio/miniconda3

LABEL author="robert mcelroy"
LABEL description="WES somatic pipeline environment"

WORKDIR /pipeline
COPY environment.yaml .

RUN conda install -n base -c conda-forge mamba && \
    mamba env update -n base -f environment.yaml && \
    conda clean -afy

ENV PATH /opt/conda/bin:$PATH
CMD ["snakemake", "--help"]
