### Filtering out multi-mappers, unmapped reads, and reads missing their pair
rule filter:
    input:
        "results/align/{sample}.bam"
    output:
        "results/filter/{sample}.bam"
    log:
        "logs/filter/{sample}.log"
    conda:
        "../envs/genomictools.yaml"
    threads: 4
    shell:
        """
        samtools view -@ {threads} -b -h -d NH:1 -q 2 -F 0x4 -F 0x8 -F 0x100 -F 0x200 -F 0x400 -F 0x800 -o {output.bam}
        """