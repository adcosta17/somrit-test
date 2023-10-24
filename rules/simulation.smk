
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

def get_maternal(wildcards):
    sample = get_sample(wildcards)
    return config["assembly_folder"]+"/"+sample+"/assembly/"+sample+".mat.fa"

def get_paternal(wildcards):
    sample = get_sample(wildcards)
    return config["assembly_folder"]+"/"+sample+"/assembly/"+sample+".pat.fa"

def get_chrom_list(wildcards):
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def get_chrom_lengths(wildcards):
    return config["chrom_lengths"]

def get_pbsim_model(wildcards):
    return config["pbsim_model"]

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



rule align_maternal:
    output:
        "HG{sample}.mat.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_maternal
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule align_paternal:
    output:
        "HG{sample}.pat.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_paternal
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule get_insert_seqs:
    input:
        mat_bam="HG{sample}.mat.sorted.bam",
        mat_bai="HG{sample}.mat.sorted.bam.bai",
        pat_bam="HG{sample}.pat.sorted.bam",
        pat_bai="HG{sample}.pat.sorted.bam.bai"
    output:
        tsv="HG{sample}/HG{sample}.insert_added_sequences.tsv", 
        mat="HG{sample}/HG{sample}.mat.assembly_subset.fa",
        pat="HG{sample}/HG{sample}.pat.assembly_subset.fa"
    threads: 1
    params:
        memory_per_thread="24G",
        repbase_seqs=get_repbase,
        mat_fasta=get_maternal,
        pat_fasta=get_paternal,
        centromeres=get_centromeres,
        chrom_lengths=get_chrom_lengths,
        script=srcdir("../scripts/get_insert_seqs.py")
    shell:
        """
        python {params.script} --input {params.repbase_seqs} --mat-bam {input.mat_bam} --pat-bam {input.pat_bam} --mat-fasta {params.mat_fasta} --pat-fasta {params.pat_fasta} --centromeres {params.centromeres} --total 500 --chrom-lengths {params.chrom_lengths} --output-folder HG{wildcards.sample} --output-prefix HG{wildcards.sample}_inserted_sequence --mat-out-fa {output.mat} --pat-out-fa {output.pat} > {output.tsv}
        """


rule simulate_seqs:
    input:
        "HG{sample}/HG{sample}.insert_added_sequences.tsv"
    output:
        tsv="HG{sample}/HG{sample}.simulated_fastq_list.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        pbsim_model=get_pbsim_model,
        script=srcdir("../scripts/generate_pbsim_run.py")
    shell:
        """
        python {params.script} --input {input} --output-folder HG{wildcards.sample} --output-prefix HG{wildcards.sample}_inserted_sequence --output-tsv {output.tsv} --pbsim-model {params.pbsim_model} --pbsim-path {config[pbsim]} > run_pbsim_{wildcards.sample}.sh
        chmod +x run_pbsim_{wildcards.sample}.sh
        ./run_pbsim_{wildcards.sample}.sh
        """

rule find_spanning_reads:
    input:
        tsv="HG{sample}/HG{sample}.simulated_fastq_list.tsv",
        inserts="HG{sample}/HG{sample}.insert_added_sequences.tsv"
    output:
        fastq="HG{sample}/HG{sample}.all_spanning_reads.fastq",
        tsv="HG{sample}/HG{sample}.spanning_reads_list.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_spanning_reads.py")
    shell:
        """
        python {params.script} --input {input.tsv} --inserts {input.inserts} --fastq {output.fastq} --tsv {output.tsv}
        """

rule simulate_fastq:
    input: 
        mat="HG{sample}/HG{sample}.mat.assembly_subset.fa",
        pat="HG{sample}/HG{sample}.pat.assembly_subset.fa"
    output:
        mat="HG{sample}/HG{sample}.mat_full.fastq",
        pat="HG{sample}/HG{sample}.pat_full.fastq"
    threads: 1
    params:
        memory_per_thread="64G",
        pbsim_model=get_pbsim_model
    shell:
        """
        {config[pbsim]} {input.mat} --prefix HG{wildcards.sample}.mat_full.pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 20 --hmm_model {params.pbsim_model} --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98
        rm HG{wildcards.sample}.mat_full.pbsim*.maf
        rm HG{wildcards.sample}.mat_full.pbsim*.ref
        cat HG{wildcards.sample}.mat_full.pbsim* >> HG{wildcards.sample}/HG{wildcards.sample}.mat_full.fastq
        rm HG{wildcards.sample}.mat_full.pbsim*
        {config[pbsim]} {input.pat} --prefix HG{wildcards.sample}.pat_full.pbsim --difference-ratio 23:31:46 --seed 10 --length-min 1000 --depth 20 --hmm_model {params.pbsim_model} --length-mean 30000 --accuracy-mean 0.95 --accuracy-min 0.85 --accuracy-max 0.98
        rm HG{wildcards.sample}.pat_full.pbsim*.maf
        rm HG{wildcards.sample}.pat_full.pbsim*.ref
        cat HG{wildcards.sample}.pat_full.pbsim* >> HG{wildcards.sample}/HG{wildcards.sample}.pat_full.fastq
        rm HG{wildcards.sample}.pat_full.pbsim*
        """

rule zip_fastq:
    input:
        "{prefix}.fastq"
    output:
        temp("{prefix}.fastq.gz")
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "bgzip -i {input}"


rule spike_fastqs:
    input:
        mat="HG{sample}/HG{sample}.mat_full.fastq.gz",
        pat="HG{sample}/HG{sample}.pat_full.fastq.gz",
        fastq="HG{sample}/HG{sample}.all_spanning_reads.fastq",
        tsv="HG{sample}/HG{sample}.spanning_reads_list.tsv"
    output:
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq",
        tsv="HG{sample}/HG{sample}.spike_in.{rep}.insert_list.txt"
    params:
        memory_per_thread="36G",
        script=srcdir("../scripts/spike_in_fastqs.py")
    shell:
        """
        python {params.script} --tsv {input.tsv} --inserts {input.fastq} --mat {input.mat} --pat {input.pat} --output-fastq {output.fastq} --rep {wildcards.rep} > {output.tsv}
        """


rule map_fastqs:
    input:
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq.gz"
    output:
        bam="HG{sample}/HG{sample}.spike_in.{rep}.bam"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.ref} {input.fastq} | samtools sort > {output.bam}
        """

rule xtea_spikein:
    input:
        bam="HG{sample}/HG{sample}.spike_in.{rep}.bam",
        bai="HG{sample}/HG{sample}.spike_in.{rep}.bam.bai"
    output:
        tsv="HG{sample}/HG{sample}.spike_in.{rep}/HG{sample}.spike_in.{rep}/classified_results.txt.Combined.txt"
    threads: 10
    benchmark:
        "benchmark/xtea/HG{sample}.spike_in.{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        cd HG{wildcards.sample}
        echo HG{wildcards.sample}.spike_in.{wildcards.rep} > HG{wildcards.sample}.spike_in.{wildcards.rep}_sample_id.txt
        echo HG{wildcards.sample}.spike_in.{wildcards.rep} ../{input.bam} > HG{wildcards.sample}.spike_in.{wildcards.rep}_long_read_bam_list.txt
        {config[xtea_dir]}/bin/xtea_long -i HG{wildcards.sample}.spike_in.{wildcards.rep}_sample_id.txt -b HG{wildcards.sample}.spike_in.{wildcards.rep}_long_read_bam_list.txt -p HG{wildcards.sample}.spike_in.{wildcards.rep} -o HG{wildcards.sample}.spike_in.{wildcards.rep}_submit_jobs.sh --clean --rmsk {config[xtea_dir]}/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r {params.ref} --cns {config[xtea_dir]}/consensus/LINE1.fa --rep {config[xtea_dir]} --xtea {config[xtea_dir]}/xtea_long/ -f 31 -y 15 -n 10 -m 190
        rm HG{wildcards.sample}.spike_in.{wildcards.rep}_submit_jobs.sh
        cd HG{wildcards.sample}.spike_in.{wildcards.rep}/HG{wildcards.sample}.spike_in.{wildcards.rep}
        chmod +x run_xTEA_pipeline.sh
        cd ../../
        ./HG{wildcards.sample}.spike_in.{wildcards.rep}/HG{wildcards.sample}.spike_in.{wildcards.rep}/run_xTEA_pipeline.sh
        """


rule sniffles_spikein:
    input:
        bam="HG{sample}/HG{sample}.spike_in.{rep}.bam",
        bai="HG{sample}/HG{sample}.spike_in.{rep}.bam.bai"
    output:
        vcf="HG{sample}/Sniffles_HG{sample}.spike_in.{rep}/HG{sample}.spike_in.{rep}.vcf"
    threads: 10
    benchmark:
        "benchmark/sniffles/HG{sample}.spike_in.{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --non-germline
        """

rule tldr_spikein:
    input:
        bam="HG{sample}/HG{sample}.spike_in.{rep}.bam",
        bai="HG{sample}/HG{sample}.spike_in.{rep}.bam.bai"
    output:
        tsv="HG{sample}/tldr_HG{sample}.spike_in.{rep}/HG{sample}.spike_in.{rep}.table.txt"
    threads: 10
    benchmark:
        "benchmark/tldr/HG{sample}.spike_in.{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        cd HG{wildcards.sample}/tldr_HG{wildcards.sample}.spike_in.{wildcards.rep}
        {config[tldr_dir]}/tldr/tldr -b ../../{input.bam} -r {params.ref} -e {config[tldr_dir]}/ref/teref.ont.human.fa -m 1 -p {threads} --flanksize 1000 --min_te_len 100
        """


rule extract_inserts_spikein:
    input:
        bam="HG{sample}/HG{sample}.spike_in.{rep}.bam",
        bai="HG{sample}/HG{sample}.spike_in.{rep}.bam.bai",
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq.gz"
    output:
        tsv="HG{sample}/HG{sample}.spike_in.{rep}.tsv",
        merged=temp("HG{sample}/HG{sample}.spike_in.{rep}.merged.txt")
    threads: 10
    benchmark:
        "benchmark/extract/HG{sample}.spike_in.{rep}.txt"
    params:
        memory_per_thread="8G",
        base_dir=get_base_dir
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-merged {params.base_dir}/{output.merged} --output-tsv {params.base_dir}/{output.tsv} --fastq-file {params.base_dir}/{input.fastq} --threads {threads} --min-flank-size 1000 --min-read-len 2000
        cd {params.base_dir}
        """

rule realign_inserts_spikein:
    input:
        bam="HG{sample}/HG{sample}.spike_in.{rep}.bam",
        tsv="HG{sample}/HG{sample}.spike_in.{rep}.tsv",
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq.gz"
    output:
        bam="HG{sample}/{chrom}.Realigned_HG{sample}.spike_in.{rep}.bam",
        tsv="HG{sample}/{chrom}.Realigned_HG{sample}.spike_in.{rep}.tsv"
    threads: 10
    benchmark:
        "benchmark/realign/HG{sample}.{chrom}.spike_in.{rep}.txt"
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        depth=max_depth,
        ref=get_ref
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.base_dir}/{input.bam} --tsv-list {params.base_dir}/{input.tsv} --fastq-list {params.base_dir}/{input.fastq} --output-dir {params.base_dir}/HG{wildcards.sample} --tsv-prefix Realigned_HG{wildcards.sample}.spike_in.{wildcards.rep} --bam-prefix Realigned_HG{wildcards.sample}.spike_in.{wildcards.rep} --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 10000 --max-depth {params.depth} --chromosome-list {wildcards.chrom} --only-realign
        cd {params.base_dir}
        """


rule merge_realign_inserts_spikein:
    input:
        bam="HG{sample}/HG{sample}.spike_in.{rep}.bam",
        cbams = expand("HG{{sample}}/{chrom}.Realigned_HG{{sample}}.spike_in.{{rep}}.bam",chrom=config["chroms"]),
        tsv = expand("HG{{sample}}/{chrom}.Realigned_HG{{sample}}.spike_in.{{rep}}.tsv",chrom=config["chroms"])
    output:
        bam="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.bam",
        tsv="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.tsv"
    threads: 1
    benchmark:
        "benchmark/merge/HG{sample}.spike_in.{rep}.txt"
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        depth=max_depth,
        ref=get_ref
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py merge --output-dir {params.base_dir}/HG{wildcards.sample} --bam-list {params.base_dir}/{input.bam} --tsv-prefix Realigned_HG{wildcards.sample}.spike_in.{wildcards.rep} --bam-prefix Realigned_HG{wildcards.sample}.spike_in.{wildcards.rep} 
        cd {params.base_dir}
        """

rule classify_inserts_spikein:
    input:
        tsv="HG{sample}/HG{sample}.spike_in.{rep}.tsv",
        bam="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.bam",
        bai="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.bam.bai",
        realign_tsv="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.tsv",
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq.gz"
    output:
        tsv=temp("HG{sample}/Realigned_classified_HG{sample}.spike_in.{rep}.tsv")
    threads: 1
    benchmark:
        "benchmark/classify/HG{sample}.spike_in.{rep}.txt"
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        repbase=get_repbase
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py classify --bam-list {params.base_dir}/{input.bam} --sample-list HG{wildcards.sample}.spike_in.{wildcards.rep} --tsv-list {params.base_dir}/{input.tsv} --realign-tsv {params.base_dir}/{input.realign_tsv} --annotation-file {params.repbase} --fastq-list {params.base_dir}/{input.fastq} --output-tsv {params.base_dir}/{output.tsv}
        cd {params.base_dir}
        """

rule filter_inserts_spikein:
    input:
        tsv="HG{sample}/Realigned_classified_HG{sample}.spike_in.{rep}.tsv",
        bam="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.bam",
        bai="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.bam.bai",
        fastq="HG{sample}/HG{sample}.spike_in.{rep}.fastq.gz"
    output:
        tsv="HG{sample}/Realigned_classified_filtered_HG{sample}.spike_in.{rep}.tsv"
    threads: 10
    benchmark:
        "benchmark/filter/HG{sample}.spike_in.{rep}.txt"
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        ref=get_ref,
        repbase=get_repbase,
        centromeres=get_centromeres,
        telomeres=get_telomeres
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py filter --threads {threads} --input-tsv {params.base_dir}/{input.tsv} --bam {params.base_dir}/{input.bam} --fastq-list {params.base_dir}/{input.fastq} --reference-genome {params.ref} --centromeres {params.centromeres} --telomeres {params.telomeres} --output-tsv {params.base_dir}/{output.tsv}
        cd {params.base_dir}
        """


rule tldr_realign_spikein:
    input:
        bam="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.bam",
        bai="HG{sample}/Realigned_HG{sample}.spike_in.{rep}.bam.bai",
    output:
        tsv="HG{sample}/Realign_tldr_HG{sample}.spike_in.{rep}/Realigned_HG{sample}.spike_in.{rep}.table.txt"
    threads: 10
    benchmark:
        "benchmark/tldr/Realign_HG{sample}.spike_in.{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        cd HG{wildcards.sample}/Realign_tldr_HG{wildcards.sample}.spike_in.{wildcards.rep}
        {config[tldr_dir]}/tldr/tldr -b ../../{input.bam} -r {params.ref} -e {config[tldr_dir]}/ref/teref.ont.human.fa -m 1 -p {threads} --flanksize 1000 --min_te_len 100
        """



rule summarize_data_spikein:
    input:
        realign_tsv="HG{sample}/Realigned_classified_filtered_HG{sample}.spike_in.{rep}.tsv",
        tldr_table="HG{sample}/tldr_HG{sample}.spike_in.{rep}/HG{sample}.spike_in.{rep}.table.txt",
        realign_tldr_table="HG{sample}/Realign_tldr_HG{sample}.spike_in.{rep}/Realigned_HG{sample}.spike_in.{rep}.table.txt",
        sniffles_vcf="HG{sample}/Sniffles_HG{sample}.spike_in.{rep}/HG{sample}.spike_in.{rep}.vcf",
        xtea_table="HG{sample}/HG{sample}.spike_in.{rep}/HG{sample}.spike_in.{rep}/classified_results.txt.Combined.txt",
        inserts_tsv="HG{sample}/HG{sample}.insert_added_sequences.tsv",
        spike_in_reads="HG{sample}/HG{sample}.spike_in.{rep}.insert_list.txt"
    output:
        tsv="simulation_results_HG{sample}.spike_in.{rep}.txt",
        tsv_500="simulation_results_500_HG{sample}.spike_in.{rep}.txt",
        tsv_2000="simulation_results_2000_HG{sample}.spike_in.{rep}.txt",
        tsv_6000="simulation_results_6000_HG{sample}.spike_in.{rep}.txt",
        tsv_LINE="simulation_results_LINE_HG{sample}.spike_in.{rep}.txt",
        tsv_Alu="simulation_results_Alu_HG{sample}.spike_in.{rep}.txt",
        tsv_SVA="simulation_results_SVA_HG{sample}.spike_in.{rep}.txt",
        benchmark="simulation_benchmark_HG{sample}.spike_in.{rep}.txt"
    threads: 1
    params:
        memory_per_thread="64G",
        somrit_test_script = srcdir("../scripts/get_simulation_metrics.py"),
        benchmark_script = srcdir("../scripts/summarize_simulation_benchmark.py"),
        xtea_folder="HG{sample}/HG{sample}.spike_in.{rep}/HG{sample}.spike_in.{rep}/classified_results.txt"
    shell:
        """
        python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {params.xtea_folder} --tldr-realign {input.realign_tldr_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} > {output.tsv}
        python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {params.xtea_folder} --tldr-realign {input.realign_tldr_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --max-size 500 > {output.tsv_500}
        python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {params.xtea_folder} --tldr-realign {input.realign_tldr_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --min-size 500 --max-size 2000 > {output.tsv_2000}
        python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {params.xtea_folder} --tldr-realign {input.realign_tldr_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --min-size 2000 > {output.tsv_6000}
        python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {params.xtea_folder} --tldr-realign {input.realign_tldr_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --family LINE > {output.tsv_LINE}
        python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {params.xtea_folder} --tldr-realign {input.realign_tldr_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --family Alu > {output.tsv_Alu}
        python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {params.xtea_folder} --tldr-realign {input.realign_tldr_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --family SVA > {output.tsv_SVA}
        python {params.benchmark_script} --rep {wildcards.rep} --sample {wildcards.sample} > {output.benchmark}
        """
