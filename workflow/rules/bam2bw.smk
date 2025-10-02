### For paired-end experiments, keep only the read that "matters" 
# Read 2 5'-end is TSS location in PRO-cap (typically)
# Read 1 5'-end is TSS location in RAMPAGE (typically)
rule get_informative_read:
    input:
        "results/filter/{sample}.bam",
    output:
        "results/filter/{sample}_informative.bam",
    log:
        "logs/get_informative_read/{sample}.log",
    conda:
        "../envs/genomictools.yaml",
    params:
        method = lookup(
            within=samples,
            query="sample_name == '{sample}'"
            cols="method",
            default="procap"
        ),
    shell:
        """
        method="{params.method}"

        if [[ "$method" == "procap" ]]; then
            samtools view -hb -f 128 {input} -o {output}

        elif [["$method" == "rampage"]]; then
            samtools view -hb -f 64 {input} -o {output}
        """

### Convert to TSS location tracks
rule bam2bg:
    input:
        get_informative_bam(),
    output:
        temp("results/bam2bg/{sample}_{strand}.bg"),
    log:
        "logs/bam2bg/{sample}_{strand}.log",
    conda:
        "../envs/genomictools.yaml",
    params:
        method = lookup(
            within=samples,
            query="sample_name == '{sample}'"
            cols="method",
            default="procap"
        ),
    shell:
        """
        method="{params.method}"
        strand="{wildcards.strand}"

        if [[ "$strand" == "plus" ]]; then
            strand_symbol="+"
        elif [[ "$strand" == "minus" ]]; then
            strand_symbol="-"

        if [[ "$method" == "procap" ]]; then
            genomeCoverageBed -ibam "$processed_bam" -bg -strand "$strand_symbol" -5 | LC_COLLATE=C sort -k1,1 -k2,2n > "$tmp_bg"
        elif [[ "$method" == "cage" ]]; then
            genomeCoverageBed -ibam "$processed_bam" -bg -strand "$strand_symbol" -5 | LC_COLLATE=C sort -k1,1 -k2,2n > "$tmp_bg"
        elif [[ "$method" == "rampage" ]]; then
            genomeCoverageBed -ibam "$processed_bam" -bg -strand "$strand_symbol" -5 | LC_COLLATE=C sort -k1,1 -k2,2n > "$tmp_bg"
        """

### Convert to bigWig
rule bg2bw:
    input:
        bg="results/bam2bg/{sample}_{strand}.bg",
        chrom=config.get("chrom_sizes"),
    output:
        "results/bg2bw/{sample}_{strand}.bigWig"
    conda:
        "../envs/genomictools.yaml"
    shell:
        """
        bedGraphToBigWig {input.bg} {input.chrom} {output}
        """


