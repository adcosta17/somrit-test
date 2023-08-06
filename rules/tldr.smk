#
# Data & Fastq related rules
#

def get_ref(wildcards):
    return config["reference"]

def get_base_dir(wildcards):
    return config["base_dir"]

def get_fastq(wildcards):
    return config["fastq"]

rule tldr_before:
    input:
        bam="HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_{frac}_{rep}.bam.bai"
    output:
        tsv="tldr_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.table.txt"
    threads: 10
    benchmark:
        "benchmark/tldr/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        cd tldr_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}
        {config[tldr_dir]}/tldr/tldr -b ../{input.bam} -r {params.ref} -e {config[tldr_dir]}/ref/teref.ont.human.fa -p {threads} --min_te_len 100
        """

rule tldr_after:
    input:
        bam="Realigned_HG{sample}_{frac}_{rep}.bam",
        bai="Realigned_HG{sample}_{frac}_{rep}.bam.bai"
    output:
        tsv="tldr_Realigned_HG{sample}_{frac}_{rep}/Realigned_HG{sample}_{frac}_{rep}.table.txt"
    threads: 10
    benchmark:
        "benchmark/tldr/Realigned_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        cd tldr_Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}
        {config[tldr_dir]}/tldr/tldr -b ../{input.bam} -r {params.ref} -e {config[tldr_dir]}/ref/teref.ont.human.fa -p {threads} --min_te_len 100
        """

