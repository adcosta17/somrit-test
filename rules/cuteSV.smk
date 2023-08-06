#
# Data & Fastq related rules
#

def get_ref(wildcards):
    return config["reference"]

def get_base_dir(wildcards):
    return config["base_dir"]

def get_fastq(wildcards):
    return config["fastq"]

rule cutesv_before:
    input:
        bam="HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_{frac}_{rep}.bam.bai"
    output:
        vcf="CuteSV_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.vcf"
    threads: 10
    log:
        "logs/cutesv/HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/cutesv/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        cuteSV {input.bam} {params.ref} {output.vcf} CuteSV_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --threads {threads} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --report_readid &> {log}
        """

rule cutesv_after:
    input:
        bam="Realigned_HG{sample}_{frac}_{rep}.bam",
        bai="Realigned_HG{sample}_{frac}_{rep}.bam.bai",
    output:
        vcf="CuteSV_Realigned_HG{sample}_{frac}_{rep}/Realigned_HG{sample}_{frac}_{rep}.vcf"
    threads: 5
    log:
        "logs/cutesv/Realigned_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/cutesv/Realigned_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="10G",
        ref=get_ref
    shell:
        """
        cuteSV {input.bam} {params.ref} {output.vcf} CuteSV_Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --threads {threads} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --report_readid &> {log}
        """