#
# Mapping related rules
#


rule sort_sams:
    output:
        bam="HG002_giab.realign.{chrom}.bam"
    threads: 5
    params:
        memory_per_thread="5G",
        sam="HG002_giab.realign.tmp.sam.{chrom}"
    shell:
        """
        samtools sort {params.sam} --threads {threads} > {output.bam}
        """

rule merge_bams:
    input:
        bams=expand("HG002_giab.realign.{c}.bam", c=config["chroms"])
    output:
        bam="HG002_giab.sorted.realigned.bam"
    threads: 5
    params:
        memory_per_thread="10G"
    shell:
        """
        samtools merge --threads {threads} {output.bam} HG002_giab.realign.*.bam
        """

## Index a bam
rule make_bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "samtools index {input}"

