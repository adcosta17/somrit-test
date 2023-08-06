
# Delete sequence to simulate somatic insertions

def get_sample(wildcards):
    return "HG00"+wildcards.sample

def get_fastq(wildcards):
    return config["fastq"]

def get_ref(wildcards):
    return config["reference"]

def max_depth(wildcards):
    return 500

def get_base_dir(wildcards):
    return config["base_dir"]

def get_repbase(wildcards):
    return config["repbase"]

def get_centromeres(wildcards):
    return config["centromeres"]

def get_telomeres(wildcards):
    return config["telomeres"]

def get_chrom_list(wildcards):
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def get_pbsim_model(wildcards):
    return config["pbsim_model"]

def get_repeat_master_grch38(wildcards):
    return config["repeat_master"]

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



rule simulate_seqs:
    output:
        tsv="simulated_fastq_list.tsv",
        data="simulated_fastq_data.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        pbsim_model=get_pbsim_model,
        script=srcdir("../scripts/generate_pbsim_run_nested.py"),
        ref=get_ref,
        rpmk=get_repeat_master_grch38,
        centromeres=get_centromeres,
        telomeres=get_telomeres,
        repbase=get_repbase
    shell:
        """
        python {params.script} --ref {params.ref} --input {params.rpmk} --output-folder simulated_fastqs --output-prefix nested --output-tsv {output.tsv} --pbsim-model {params.pbsim_model} --pbsim-path {config[pbsim]} --output-data {output.data} --centromeres {params.centromeres} --telomeres {params.telomeres} --seqs {params.repbase} --rpmk {params.rpmk} > run_pbsim_nested.sh
        chmod +x run_pbsim_nested.sh
        ./run_pbsim_nested.sh
        """

rule find_spanning_reads:
    input:
        tsv="simulated_fastq_list.tsv"
    output:
        tsv="nested_insert_comp.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_spanning_reads_nested.py")
    shell:
        """
        python {params.script} --input {input.tsv} --tsv {output.tsv}
        """


