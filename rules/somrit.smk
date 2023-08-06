#
# Data & Fastq related rules
#

def max_depth(wildcards):
    return int(float(wildcards.frac)*600)

def get_ref(wildcards):
    return config["reference"]

def get_base_dir(wildcards):
    return config["base_dir"]

def get_fastq(wildcards):
    return config["fastq"]

def get_repbase(wildcards):
    return config["repbase"]

def get_centromeres(wildcards):
    return config["centromeres"]

def get_telomeres(wildcards):
    return config["telomeres"]

rule extract_inserts:
    input:
        bam="HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_{frac}_{rep}.bam.bai",
        fastq="HG{sample}_{frac}_{rep}.fastq.bgz"
    output:
        tsv="HG{sample}_{frac}_{rep}.tsv",
        merged=temp("HG{sample}_merged_{frac}_{rep}.txt")
    threads: 10
    log:
        "logs/somrit/HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="8G",
        base_dir=get_base_dir
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-merged {params.base_dir}/{output.merged} --output-tsv {params.base_dir}/{output.tsv} --fastq-file {params.base_dir}/{input.fastq} --threads {threads} &> {params.base_dir}/{log}
        cd {params.base_dir}
        """

rule realign_inserts:
    input:
        tsv="HG{sample}_{frac}_{rep}.tsv",
        fastq="HG{sample}_{frac}_{rep}.fastq.bgz"
    output:
        bam="{chrom}.Realigned_HG{sample}_{frac}_{rep}.bam",
        tsv="{chrom}.Realigned_HG{sample}_{frac}_{rep}.tsv"
    threads: 10
    log:
        "logs/somrit/{chrom}.Realigned_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit/{chrom}.Realigned_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        depth=max_depth,
        ref=get_ref,
        bam="HG{sample}_{frac}_{rep}.bam",
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.base_dir}/{params.bam} --tsv-list {params.base_dir}/{input.tsv} --fastq-list {params.base_dir}/{input.fastq} --output-dir {params.base_dir} --tsv-prefix Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --bam-prefix Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 10000 --max-depth {params.depth} --chromosome-list {wildcards.chrom} --only-realign &> {params.base_dir}/{log}
        cd {params.base_dir}
        """


rule merge_realign_inserts:
    input:
        cbams = expand("{chrom}.Realigned_HG{{sample}}_{{frac}}_{{rep}}.bam",chrom=config["chroms"]),
        tsv = expand("{chrom}.Realigned_HG{{sample}}_{{frac}}_{{rep}}.tsv",chrom=config["chroms"])
    output:
        bam="Realigned_HG{sample}_{frac}_{rep}.bam",
        tsv="Realigned_HG{sample}_{frac}_{rep}.tsv"
    threads: 1
    log:
        "logs/somrit/Realigned_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit/Realigned_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        depth=max_depth,
        ref=get_ref,
        bam="HG{sample}_{frac}_{rep}.bam",
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py merge --output-dir {params.base_dir} --bam-list {params.base_dir}/{params.bam} --tsv-prefix Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --bam-prefix Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} &> {params.base_dir}/{log}
        cd {params.base_dir}
        """

rule classify_inserts:
    input:
        realign_tsv="Realigned_HG{sample}_{frac}_{rep}.tsv",
    output:
        tsv="Realigned_classified_HG{sample}_{frac}_{rep}.tsv"
    threads: 1
    log:
        "logs/somrit/Realigned_classified_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit/Realigned_classified_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        tsv="HG{sample}_{frac}_{rep}.tsv",
        bam="Realigned_HG{sample}_{frac}_{rep}.bam",
        bai="Realigned_HG{sample}_{frac}_{rep}.bam.bai",
        fastq="HG{sample}_{frac}_{rep}.fastq.bgz",
        repbase=get_repbase
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py classify --bam-list {params.base_dir}/{params.bam} --sample-list HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --tsv-list {params.base_dir}/{params.tsv} --realign-tsv {params.base_dir}/{input.realign_tsv} --annotation-file {params.repbase} --fastq-list {params.base_dir}/{params.fastq} --output-tsv {params.base_dir}/{output.tsv} &> {params.base_dir}/{log}
        cd {params.base_dir}
        """

rule filter_inserts:
    input:
        tsv="Realigned_classified_HG{sample}_{frac}_{rep}.tsv"
    output:
        tsv="Realigned_classified_filtered_HG{sample}_{frac}_{rep}.tsv"
    threads: 10
    log:
        "logs/somrit/Realigned_classified_filtered_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit/Realigned_classified_filtered_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        ref=get_ref,
        repbase=get_repbase,
        centromeres=get_centromeres,
        telomeres=get_telomeres,
        bam="Realigned_HG{sample}_{frac}_{rep}.bam",
        bai="Realigned_HG{sample}_{frac}_{rep}.bam.bai"
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py filter --threads {threads} --input-tsv {params.base_dir}/{input.tsv} --bam {params.base_dir}/{params.bam} --fastq-list NA --reference-genome {params.ref} --centromeres {params.centromeres} --telomeres {params.telomeres} --output-tsv {params.base_dir}/{output.tsv} --min-reads 3 &> {params.base_dir}/{log}
        cd {params.base_dir}
        """



