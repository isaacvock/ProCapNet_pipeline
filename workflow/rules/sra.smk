rule get_sra_pe:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "v7.2.0/bio/sra-tools/fasterq-dump"



rule get_sra_se:
    output:
        "sra/{accession}.fastq",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "v7.2.0/bio/sra-tools/fasterq-dump"