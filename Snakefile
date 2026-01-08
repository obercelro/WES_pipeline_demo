# WES analysis pipeline
# author: robert mcelroy

configfile: "config/config.yaml"

import pandas as pd
import os
import yaml

# ----------
# dynamic resource loading
# allows switching between laptop (low RAM) and HPC (high RAM) without changing code

resource_config = config.get("resources", "config/resources.yaml")
if os.path.exists(resource_config):
    with open(resource_config, 'r') as f:
        RES = yaml.safe_load(f)["rules"]
else:
    # fallback defaults
    RES = {
        "bwa_mem": {"threads": 4, "mem_mb": 8192, "time": "08:00:00"},
        "mark_duplicates": {"mem_mb": 4096},
        "bqsr": {"mem_mb": 4096},
        "mutect2": {"mem_mb": 8192},
        "haplotype_caller": {"mem_mb": 8192},
        "manta_sv": {"threads": 4, "mem_mb": 4096}
    }

# ----------
# helper functions

# if keep_intermediates is False in config, files are marked as temp()
# usage: output: intermediate("path/to/file")
def intermediate(path):
    if config.get("keep_intermediates", True):
        return path
    return temp(path)

# parse input csvs
units = pd.read_csv(config["units"], dtype=str).set_index(["sample_id", "lane"], drop=False)
samples = pd.read_csv(config["samples"], dtype=str).set_index("sample_id", drop=False)

def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.lane), ["fastq_1", "fastq_2"]].dropna()

def get_normal_bam(wildcards):
    # lookup the normal sample ID associated with this patient
    patient = samples.loc[wildcards.sample, "patient_id"]
    normal_id = samples[(samples["patient_id"] == patient) & (samples["type"] == "normal")]["sample_id"].item()
    return f"results/recal/{normal_id}.recal.bam"

# ----------
# target outputs

rule all:
    input:
        "results/qc/multiqc_report.html",
        # somatic variants
        expand("results/vcf/somatic/annotated/{sample}.somatic.ann.vcf.gz", sample=samples[samples["type"] == "tumor"]["sample_id"]),
        # germline variants
        expand("results/vcf/germline/{sample}.germline.vcf.gz", sample=samples[samples["type"] == "normal"]["sample_id"]),
        # structural variants
        expand("results/sv/manta/{sample}/results/variants/somaticSV.vcf.gz", sample=samples[samples["type"] == "tumor"]["sample_id"]),
        # tabular analysis files
        expand("results/maf/{sample}.maf", sample=samples[samples["type"] == "tumor"]["sample_id"])
        # NOTE: cram compression is optional and not listed here. 
        # run 'snakemake --use-conda results/archived/{sample}.cram' to trigger it.

# ----------
# upstream processing

rule fastp:
    input: get_fastq
    output:
        r1 = "results/trimmed/{sample}_{lane}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_{lane}_R2.fastq.gz",
        json = "results/qc/fastp/{sample}_{lane}.json",
        html = "results/qc/fastp/{sample}_{lane}.html"
    resources: mem_mb=4000
    threads: 4
    shell:
        "fastp -i {input[0]} -I {input[1]} -o {output.r1} -O {output.r2} -j {output.json} -h {output.html} -w {threads}"

rule bwa_index:
    input:
        config["ref"]["genome"]
    output:
        multiext(config["ref"]["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    resources:
        mem_mb=16000
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        r1 = "results/trimmed/{sample}_{lane}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_{lane}_R2.fastq.gz",
        ref = config["ref"]["genome"],
        idx = multiext(config["ref"]["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output: 
        intermediate("results/mapped/{sample}_{lane}.sorted.bam")
    threads: RES["bwa_mem"]["threads"]
    resources: 
        mem_mb = RES["bwa_mem"]["mem_mb"],
        time = RES["bwa_mem"]["time"]
    shell:
        "bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}_{wildcards.lane}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' "
        "{input.ref} {input.r1} {input.r2} | samtools sort -@ {threads} -m 1G -o {output} -"

rule mark_duplicates:
    input: lambda wc: expand("results/mapped/{sample}_{lane}.sorted.bam", sample=wc.sample, lane=units.loc[wc.sample, "lane"])
    output:
        bam = intermediate("results/dedup/{sample}.dedup.bam"),
        metrics = "results/qc/dedup/{sample}.metrics.txt"
    resources: mem_mb=RES["mark_duplicates"]["mem_mb"]
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} --CREATE_INDEX true"

rule bqsr:
    input:
        bam = "results/dedup/{sample}.dedup.bam",
        ref = config["ref"]["genome"],
        known = config["ref"]["known_variants"]
    output:
        table = "results/recal/{sample}.recal_data.table",
        bam = "results/recal/{sample}.recal.bam"
    resources: mem_mb=RES["bqsr"]["mem_mb"]
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known} -O {output.table}
        gatk ApplyBQSR -I {input.bam} -R {input.ref} --bqsr-recal-file {output.table} -O {output.bam}
        """

# ----------
# somatic variant calling

rule mutect2:
    input:
        tumor = "results/recal/{sample}.recal.bam",
        normal = get_normal_bam,
        ref = config["ref"]["genome"],
        intervals = config["ref"]["intervals"]
    output:
        vcf = "results/vcf/somatic/raw/{sample}.somatic.vcf.gz",
        f1r2 = "results/vcf/somatic/raw/{sample}.f1r2.tar.gz"
    resources: mem_mb=RES["mutect2"]["mem_mb"]
    shell:
        """
        gatk Mutect2 -R {input.ref} -I {input.tumor} -I {input.normal} \
            -normal $(basename {input.normal} .recal.bam) \
            -L {input.intervals} --f1r2-tar-gz {output.f1r2} -O {output.vcf}
        """

rule filter_mutect:
    input:
        vcf = "results/vcf/somatic/raw/{sample}.somatic.vcf.gz",
        ref = config["ref"]["genome"]
    output: "results/vcf/somatic/filtered/{sample}.somatic.filtered.vcf.gz"
    shell: "gatk FilterMutectCalls -R {input.ref} -V {input.vcf} -O {output}"

rule vep_annotate:
    input:
        vcf = "results/vcf/somatic/filtered/{sample}.somatic.filtered.vcf.gz",
        cache = config["ref"]["vep_cache"]
    output: "results/vcf/somatic/annotated/{sample}.somatic.ann.vcf.gz"
    shell:
        "vep --cache --dir_cache {input.cache} --assembly GRCh38 --offline --input_file {input.vcf} --output_file {output} --format vcf --vcf --compress_output gzip --force_overwrite"

# ----------
# germline & sv calling

rule haplotype_caller:
    input:
        bam = "results/recal/{sample}.recal.bam",
        ref = config["ref"]["genome"],
        intervals = config["ref"]["intervals"]
    output: "results/vcf/germline/{sample}.germline.vcf.gz"
    resources: mem_mb=RES["haplotype_caller"]["mem_mb"]
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} -L {input.intervals} -ERC GVCF"

rule manta_sv:
    input:
        tumor = "results/recal/{sample}.recal.bam",
        normal = get_normal_bam,
        ref = config["ref"]["genome"]
    output:
        vcf = "results/sv/manta/{sample}/results/variants/somaticSV.vcf.gz"
    params:
        run_dir = "results/sv/manta/{sample}"
    conda: "envs/manta.yaml"
    threads: RES["manta_sv"]["threads"]
    resources: mem_mb=RES["manta_sv"]["mem_mb"]
    shell:
        """
        rm -rf {params.run_dir}
        configManta.py --normalBam {input.normal} --tumorBam {input.tumor} --referenceFasta {input.ref} --runDir {params.run_dir} --exome
        {params.run_dir}/runWorkflow.py -m local -j {threads}
        """

# ---------- 
# analysis & reporting

rule vcf2maf:
    input:
        vcf = "results/vcf/somatic/annotated/{sample}.somatic.ann.vcf.gz",
        ref = config["ref"]["genome"]
    output:
        maf = "results/maf/{sample}.maf"
    params:
        tumor_id = "{sample}",
        normal_id = lambda wc: "Normal_1000G"
    conda: "envs/maf.yaml"
    shell:
        """
        gunzip -c {input.vcf} > {input.vcf}.temp
        
        vcf2maf.pl \
            --input-vcf {input.vcf}.temp \
            --output-maf {output.maf} \
            --tumor-id {params.tumor_id} \
            --normal-id {params.normal_id} \
            --ref-fasta {input.ref} \
            --ncbi-build GRCh38 \
            --inhibit-vep
        
        rm {input.vcf}.temp
        """

rule samtools_stats:
    input: "results/dedup/{sample}.dedup.bam"
    output: "results/qc/bam/{sample}.stats"
    shell: "samtools stats {input} > {output}"

rule multiqc:
    input:
        expand("results/qc/fastp/{u.sample_id}_{u.lane}.json", u=units.itertuples()),
        expand("results/qc/dedup/{sample}.metrics.txt", sample=samples["sample_id"]),
        expand("results/qc/bam/{sample}.stats", sample=samples["sample_id"])
    output: "results/qc/multiqc_report.html"
    shell: "multiqc results/qc -o results/qc -f"

# ----------
# optional: archiving (cram)

rule bam_to_cram:
    input:
        bam = "results/recal/{sample}.recal.bam",
        ref = config["ref"]["genome"]
    output:
        cram = "results/archived/{sample}.cram",
        crai = "results/archived/{sample}.cram.crai"
    threads: 4
    shell:
        """
        samtools view -C -T {input.ref} -o {output.cram} {input.bam}
        samtools index {output.cram}
        """
