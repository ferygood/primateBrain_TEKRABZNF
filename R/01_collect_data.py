import pandas as pd

# read sra files, accession PRJNA527986
df = pd.read_csv("SraRunTalbe.txt", sep=",", header=0)

# append all accession number to a list
SAMPLES = []

for i in df['Run']:
    SAMPLES.append(i)

#############################    
# snakemake workflows #######
#############################

# 1. fasterq-dump.snake
# getting fastq files for Primate Brain Data
rule all:
    input:
        expand("<path>/fastq/{accession}_1.fastq.gz", accession=SAMPLES)
        
rule get_fastq_pe_gz:
    output:
        # the wildcard name must be accession
        "<path>/fastq/{accession}_1.fastq.gz",
        "<path>/fastq/{accession}_2.fastq.gz",
    log:
        "<path>/logs/{accession}.log"
    params:
        extra="--skip-technical"
    threads: 20
    wrapper:
        "v1.7.2/bio/sra-tools/fasterq-dump"
        
# 2. Download fastq expression from Synapse
# getting fastq files for Mayo Data
"""
# in bash command
synapse -u <account> -p <token> get -r syn8612203 --downloadLocation ./
synapse -u <account> -p <token> get -r syn8612213 --downloadLocation ./
"""

        
        
