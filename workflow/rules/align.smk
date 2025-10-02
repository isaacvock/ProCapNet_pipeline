if ALIGNER == "star":

    # Build STAR index
    rule index:
        input:
            fasta=config["genome"],
            gtf=config["annotation"],
        output:
            outdir=directory(config.get("indices", "star_index")),
            outfile=star_indices,
        params:
            extra=config["star_index_params"],
        log:
            "logs/index/star_index_genome.log",
        conda:
            "../envs/star.yaml"
        threads: 24
        script:
            "../scripts/alignment/star-index.py"

    # Align with STAR
        # I think passing GTF here is redundant given how
        # indexing is done, but I don't think it hurts and 
        # can catch provided-index-edge-cases
    rule align:
        input:
            unpack(get_fq),
            idx=config.get("indices", "star_index"),
            gtf=config["annotation"],
        output:
            aln="results/align/{sample}.bam",
            reads_per_gene="results/align/{sample}_ReadsPerGene.out.tab",
            log="results/align/{sample}_Log.out",
            sj="results/align/{sample}_SJ.out.tab",
            log_final="results/align/{sample}_Log.final.out",
        log:
            "logs/align/{sample}.log",
        params:
            extra=lambda wc, input: " ".join(
                [
                    "--outSAMtype BAM SortedByCoordinate",
                    "--quantMode GeneCounts",
                    f'--sjdbGTFfile "{input.gtf}"',
                    lookup(within=config, dpath="star_align_extra", default=""),
                ]
            ),
        threads: 24
        wrapper:
            "v7.2.0/bio/star/align"

elif ALIGNER == "bowtie2":
    ValueError("aligner currently can only be star")

elif ALIGNER == "bwa-mem2":
    ValueError("aligner currently can only be star")

else:
    ValueError("aligner must be star, bowtie2, or bwa-mem2!!")