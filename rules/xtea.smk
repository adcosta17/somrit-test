#
# Data & Fastq related rules
#

def get_ref(wildcards):
    return config["reference"]

def get_base_dir(wildcards):
    return config["base_dir"]

def get_fastq(wildcards):
    return config["fastq"]

rule xtea_before:
    output:
        tsv="HG{sample}_{frac}_{rep}/classified_results.txt.Combined.txt"
    threads: 10
    benchmark:
        "benchmark/xtea/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref,
        bam="HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_{frac}_{rep}.bam.bai"
    shell:
        """
        echo HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} > HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_sample_id.txt
        echo HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} {params.bam} > HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_long_read_bam_list.txt
        {config[xtea_dir]}/bin/xtea_long -i HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_sample_id.txt -b HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_long_read_bam_list.txt -p . -o HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_submit_jobs.sh --clean --rmsk {config[xtea_dir]}/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r {params.ref} --cns {config[xtea_dir]}/consensus/LINE1.fa --rep {config[xtea_dir]} --xtea {config[xtea_dir]}/xtea_long/ -f 31 -y 15 -n 10 -m 190 --clean
        rm HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_submit_jobs.sh
        cd HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}
        chmod +x run_xTEA_pipeline.sh
        cd ..
        ./HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}/run_xTEA_pipeline.sh
        """

rule xtea_after:
    input:
        bam="Realigned_HG{sample}_{frac}_{rep}.bam",
        bai="Realigned_HG{sample}_{frac}_{rep}.bam.bai"
    output:
        tsv="Realigned_HG{sample}_{frac}_{rep}/classified_results.txt.Combined.txt"
    threads: 10
    benchmark:
        "benchmark/xtea/Realigned_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        echo Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} > Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_sample_id.txt
        echo Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} {input.bam} > Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_long_read_bam_list.txt
        {config[xtea_dir]}/bin/xtea_long -i Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_sample_id.txt -b Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_long_read_bam_list.txt -p . -o Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_submit_jobs.sh --clean --rmsk {config[xtea_dir]}/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r {params.ref} --cns {config[xtea_dir]}/consensus/LINE1.fa --rep {config[xtea_dir]} --xtea {config[xtea_dir]}/xtea_long/ -f 31 -y 15 -n 10 -m 190 --clean
        rm Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_submit_jobs.sh
        cd Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}
        chmod +x run_xTEA_pipeline.sh
        cd ..
        ./Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}/run_xTEA_pipeline.sh
        """

