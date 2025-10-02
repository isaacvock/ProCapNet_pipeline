### Read 2 5'-end is TSS location in PRO-cap (typically)
### Read 1 5'-end is TSS location in RAMPAGE (typically)
rule get_informative_read:
    input:
        "results/filter/{sample}.bam",
    output:
        "results/filter/{sample}_informative.bam",
    log:
        "logs/get_informative_read/{sample}.log",
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

rule bam2bg_plus:
    input:
        get_informative_bam(),
    output:
        temp("results/bam2bg/{sample}_plus.bg"),
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
            genomeCoverageBed -ibam "$processed_bam" -bg -strand "$strand_symbol" -5 | grep -v "_" | LC_COLLATE=C sort -k1,1 -k2,2n > "$tmp_bg"
        elif [[ "$method" == "cage" ]]; then
            genomeCoverageBed -ibam "$processed_bam" -bg -strand "$strand_symbol" -5 | grep -e "^chr[0-9XY]*	" | LC_COLLATE=C sort -k1,1 -k2,2n > "$tmp_bg"
        elif [[ "$method" == "rampage" ]]; then
            genomeCoverageBed -ibam "$processed_bam" -bg -strand "$strand_symbol" -5 | grep -e "^chr[0-9XY]*	" | LC_COLLATE=C sort -k1,1 -k2,2n > "$tmp_bg"
        """

# For now, I have removed Kelly's method-specific filtering, which is something
# that presumably can be done post-hoc
rule bam2bg_minus:
    input:
        get_informative_bam(),
    output:
        temp("results/bam2bg/{sample}_minus.bg"),
    conda:
        "../envs/genomictools.yaml",
    params:
        method = lookup(
            within=samples,
            query="sample_name == '{sample}'"
            cols="method",
            default="procap"
        ),
        bam2bg_filtering = config.get("bam2bg_filtering", "")
    shell:
        """
        method="{params.method}"
        if [[ "$method" == "procap" ]]; then
            genomeCoverageBed -ibam {input} -bg -strand - -5 | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        elif [[ "$method" == "cage" ]]; then
            genomeCoverageBed -ibam {input} -bg -strand - -5 | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        elif [[ "$method" == "rampage" ]]; then
            genomeCoverageBed -ibam {input} -bg -strand - -5 | LC_COLLATE=C sort -k1,1 -k2,2n > {output}
        """

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


