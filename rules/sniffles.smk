#
# Data & Fastq related rules
#

def get_ref(wildcards):
    return config["reference"]

def get_base_dir(wildcards):
    return config["base_dir"]

def get_fastq(wildcards):
    return config["fastq"]

rule sniffles_before:
    input:
        bam="HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_{frac}_{rep}.bam.bai"
    output:
        vcf="Sniffles_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.vcf"
    threads: 10
    log:
        "logs/sniffles/HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/sniffles/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        tmp="HG{sample}_{frac}_{rep}.tmp.bam",
        ref=get_ref
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --output-rnames &> {log}
        """

rule sniffles_after:
    input:
        bam="Realigned_HG{sample}_{frac}_{rep}.bam",
        bai="Realigned_HG{sample}_{frac}_{rep}.bam.bai",
    output:
        vcf="Sniffles_Realigned_HG{sample}_{frac}_{rep}/Realigned_HG{sample}_{frac}_{rep}.vcf"
    threads: 5
    log:
        "logs/sniffles/Realigned_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/sniffles/Realigned_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="10G",
        tmp="Realigned_HG{sample}_{frac}_{rep}.tmp.bam",
        ref=get_ref
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --output-rnames &> {log}
        """

