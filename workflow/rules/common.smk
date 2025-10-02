import glob

import pandas as pd
from snakemake.utils import validate


##### HELPFUL CONFIG PARSING #####

# Which adapter trimmer to use
TRIMMER = config.get("trimmer", "cutadapt")

# Need some sort of tangible file to confirm that alignment
# index has been built/exists
star_indices = glob.glob(f"{INDEX_PATH}/genomeParameters.txt")

##### FUNCTIONS USED THROUGHOUT #####

### Load sample information
samples = (
    pd.read_csv(config["samples"], dtype={"sample_name": str, "replicate_group": str})
    .set_index(["sample_name", "replicate_group"], drop=False)
    .sort_index()
)


### Get paths to fastqs for trimming
def get_units_fastqs(wildcards):
    s = samples.loc[(wildcards.sample)]
    if pd.isna(s["fq1"]):
        # SRA sample (always paired-end for now)
        accession = s["sra"]
        PE = s["PE"]
        if PE:
            return expand(
                "sra/{accession}_{read}.fastq",
                accession=accession,
                read=["R1", "R2"],
            )
        else:
            return expand(
                "sra/{accession}_{read}.fastq",
                accession=accession,
                read=["R1", "R2"],
            )
    if not is_paired_end(wildcards.sample):
        return [
            u["fq1"],
        ]
    else:
        return [u["fq1"], u["fq2"]]


### Figure out if sample is paired end
def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired



### Get paths to fastqs for alignment
def get_fq(wildcards):
    if config["trim_reads"]:
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "results/trimmed/{sample}_{read}.fastq.gz",
                        read=["R1", "R2"],
                        **wildcards,
                    ),
                )
            )
        # single end sample
        return {
            "fq1": "results/trimmed/{sample}_single.fastq.gz".format(
                **wildcards
            )
        }
    else:
        # no trimming, use raw reads
        fqs = get_units_fastqs(wildcards)
        if len(fqs) == 1:
            return {"fq1": f"{fqs[0]}"}
        elif len(fqs) == 2:
            return {"fq1": f"{fqs[0]}", "fq2": f"{fqs[1]}"}
        else:
            raise ValueError(f"Expected one or two fastq file paths, but got: {fqs}")


### Get bam files that will be used to make bedgraphs (i.e., have unique information about TSS location)
def get_informative_bam(wildcards):
    s = samples.loc[wildcards.sample]

    method = s["method"]

    if(len(method) > 1):
        ValueError("sample_names in samples CSV must be unique!")

    if method == "procap" or method == "rampage":
        return expand(
            "results/filter/{sample}_informative.bam",
            sample = wildcards.sample
        )
    else:
        return expand(
            "results/filter/{sample}.bam",
            sample = wildcards.sample
        )



### Get final output
def get_final_output():
    # What?
    # 1) bam files for CAGE-seq and DNase-seq data
    # 2) Informative bam files for PRO-cap and RAMPAGE
    # 3) +/- Bigwigs, but only for non-DNase-seq data
