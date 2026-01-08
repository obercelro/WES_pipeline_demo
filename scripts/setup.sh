#!/bin/bash
set -e

# colors
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # no color
# resource buckets
BROAD_BUCKET="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
CACHE_DIR="data/refs/vep_cache"
VEP_URL="https://ftp.ensembl.org/pub/release-106/variation/indexed_vep_cache/homo_sapiens_vep_106_GRCh38.tar.gz"

echo -e "${GREEN}  WES example data setup                                ${NC}"
echo -e "${GREEN}  Source: Broad + EBI                                   ${NC}"

mkdir -p data/refs data/fastq config results logs
mkdir -p data/refs/vep_cache

download_verified() {
    local url=$1
    local output=$2
    local file_type=$3 

    echo -n "  -> fetching $(basename $output)!"
 
    if curl -L -f -o "$output" "$url" 2>/dev/null; then
        echo -n "downloaded!"
    else
        echo -e "${RED}FAILED (HTTP Error).${NC}"
        echo "URL: $url"
        exit 1
    fi

    # make sure its not empty
    if [ ! -s "$output" ]; then
        echo -e "${RED}ERROR (empty file).${NC}"
        rm "$output"
        exit 1
    fi

    # if GZIP, check for magic number
    if [ "$file_type" == "GZ" ]; then
        if [[ $(xxd -p -l 2 "$output") != "1f8b" ]]; then
            echo -e "${RED}ERROR (corrupt / not GZIP).${NC}"
            rm "$output"
            exit 1
        fi
    fi

    echo -e "${GREEN}done!${NC}"
}

# ----------
echo "[1/4] downloading reference assets."

# genome 
download_verified "${BROAD_BUCKET}/Homo_sapiens_assembly38.fasta" "data/refs/Homo_sapiens_assembly38.fasta" "TEXT"
download_verified "${BROAD_BUCKET}/Homo_sapiens_assembly38.fasta.fai" "data/refs/Homo_sapiens_assembly38.fasta.fai" "TEXT"
download_verified "${BROAD_BUCKET}/Homo_sapiens_assembly38.dict" "data/refs/Homo_sapiens_assembly38.dict" "TEXT"

# support files
download_verified "${BROAD_BUCKET}/Homo_sapiens_assembly38.dbsnp138.vcf.gz" "data/refs/dbsnp.vcf.gz" "GZ"
download_verified "${BROAD_BUCKET}/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi" "data/refs/dbsnp.vcf.gz.tbi" "GZ"

# intervals
echo -n "  -> generating intervals (BED)."
download_verified "${BROAD_BUCKET}/wgs_calling_regions.hg38.interval_list" "data/refs/wgs_calling_regions.interval_list" "TEXT"
grep -v '@' data/refs/wgs_calling_regions.interval_list | awk '{print $1 "\t" $2 "\t" $3}' > data/refs/intervals.bed
echo -e "${GREEN}done!${NC}"

# ----------
echo "[2/4] downloading WES data."

# normal (daughter)
download_verified "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098401/SRR098401_1.fastq.gz" "data/fastq/Normal_R1.fastq.gz" "GZ"
download_verified "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098401/SRR098401_2.fastq.gz" "data/fastq/Normal_R2.fastq.gz" "GZ"

# tumor (father)
download_verified "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098359/SRR098359_1.fastq.gz" "data/fastq/Tumor_R1.fastq.gz" "GZ"
download_verified "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098359/SRR098359_2.fastq.gz" "data/fastq/Tumor_R2.fastq.gz" "GZ"

# ----------
echo "[3/4] generating configs."

cat <<EOF > config/units.csv
sample_id,lane,fastq_1,fastq_2
Normal_1000G,L001,data/fastq/Normal_R1.fastq.gz,data/fastq/Normal_R2.fastq.gz
Tumor_1000G,L001,data/fastq/Tumor_R1.fastq.gz,data/fastq/Tumor_R2.fastq.gz
EOF

cat <<EOF > config/samples.csv
patient_id,sample_id,type
PATIENT_001,Normal_1000G,normal
PATIENT_001,Tumor_1000G,tumor
EOF

# check if resources.yaml exists (otherwise a default is created)
if [ ! -f config/resources.yaml ]; then
cat <<EOF > config/resources.yaml
default:
  threads: 1
  mem_mb: 1024
  time: "01:00:00"
rules:
  bwa_mem:
    threads: 4
    mem_mb: 8192
    time: "08:00:00"
  mark_duplicates:
    mem_mb: 6144
  bqsr:
    mem_mb: 4096
  mutect2:
    threads: 2
    mem_mb: 8192
    time: "12:00:00"
  haplotype_caller:
    mem_mb: 8192
  manta_sv:
    threads: 4
    mem_mb: 4096
EOF
fi

echo "[4/4] double checking VEP cache."
(
  mkdir -p "$CACHE_DIR"
  cd "$CACHE_DIR"

  if [ ! -d "homo_sapiens" ]; then
    echo "no cache, downloading VEP 106"
    curl -f -L -O $VEP_URL
    tar -xzf homo_sapiens_vep_106_GRCh38.tar.gz
    rm homo_sapiens_vep_106_GRCh38.tar.gz
  fi
)

echo "✓✓"

echo -e "${GREEN}  setup complete! have fun!                                      ${NC}"
echo -e "${GREEN}  to run: snakemake --cores 4 --use-conda --conda-frontend conda ${NC}"
echo -e "${GREEN}  --use-conda & --conda-frontend are enabled to allow creation   ${NC}"
echo -e "${GREEN}  of separate environment for manta (snakemake doesn't like to   ${NC}"
echo -e "${GREEN}  have multiple python versions in the same environment)         ${NC}"
