
# Delete sequence to simulate somatic insertions

def get_sample(wildcards):
    return config["full_sample"]

def get_assembly(wildcards):
    sample = get_sample(wildcards)
    return config["assembly_folder"]+"/"+sample+"/assembly/"+sample+"."+wildcards.hap+".fa"

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

def get_tsv_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+"/"+sample+".tsv,"
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_sample_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += sample+","
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_fastq_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+"/fastq/"+sample+".fastq.gz,"
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_bam_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+"/mapped/"+sample+".bam,"
    tsv_list = tsv_list[:-1]
    return tsv_list

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



rule map_fastqs:
    input:
        fastq="{sample}/fastq/{sample}.fastq.gz"
    output:
        bam="{sample}/mapped/{sample}.bam"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.ref} {input.fastq} | samtools sort > {output.bam}
        """


rule tldr:
    input:
        bam=expand("{sample}/mapped/{sample}.bam", sample=config["samples"]),
        bai=expand("{sample}/mapped/{sample}.bam.bai", sample=config["samples"])
    output:
        tsv="combined/tldr_combined_{chrom}/combined.table.txt"
    threads: 10
    params:
        memory_per_thread="30G",
        ref=get_ref,
        bam_list=get_bam_list
    shell:
        """
        cd combined/tldr_combined_{wildcards.chrom}
        echo {wildcards.chrom} > chrom.txt
        {config[tldr_dir]}/tldr/tldr -b {params.bam_list} -r {params.ref} -e {config[tldr_dir]}/ref/teref.mouse.fa -n {config[tldr_dir]}/ref/nonref.collection.mm10.chr.bed.gz --keep_pickles -m 1 -p {threads} -c chrom.txt --flanksize 1000 -o combined
        """

rule merge_tldr:
    input:
        expand("combined/tldr_combined_{chrom}/combined.table.txt", chrom=config["chroms"])
    output:
        "combined/combined.tldr.table.txt"
    threads: 1
    params:
        memory_per_thread="4G",
    shell:
        """
        head -1 combined/tldr_combined_chr1/combined.table.txt > {output}
        for file in {input}
        do
          tail -n+2 "$file" >> {output}
        done
        """


rule extract_inserts:
    input:
        bam="{sample}/mapped/{sample}.bam",
        bai="{sample}/mapped/{sample}.bam.bai",
        fastq="{sample}/fastq/{sample}.fastq.gz"
    output:
        tsv="{sample}/{sample}.tsv",
        merged=temp("{sample}/{sample}.merged.txt")
    threads: 10
    params:
        memory_per_thread="8G",
        base_dir=get_base_dir
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-merged {params.base_dir}/{output.merged} --output-tsv {params.base_dir}/{output.tsv} --fastq-file {params.base_dir}/{input.fastq} --threads {threads} --min-flank-size 1000 --min-read-len 2000
        cd {params.base_dir}
        """

rule realign_inserts:
    input:
        bam=expand("{sample}/mapped/{sample}.bam",sample=config["samples"]),
        tsv=expand("{sample}/{sample}.tsv",sample=config["samples"]),
        fastq=expand("{sample}/fastq/{sample}.fastq.gz",sample=config["samples"])
    output:
        bam="combined_realign_all/{chrom}.Realigned.bam",
        tsv="combined_realign_all/{chrom}.Realigned.tsv"
    threads: 10
    params:
        memory_per_thread="70G",
        base_dir=get_base_dir,
        depth=max_depth,
        bam_list=get_bam_list,
        tsv_list=get_tsv_list,
        fastq_list=get_fastq_list,
        ref=get_ref
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.bam_list} --tsv-list {params.tsv_list} --fastq-list {params.fastq_list} --output-dir {params.base_dir}/combined_realign_all --tsv-prefix Realigned --bam-prefix Realigned --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 7000 --max-depth {params.depth} --chromosome-list {wildcards.chrom} --filter-pass --only-realign --filter-pass
        cd {params.base_dir}
        """


rule merge_realign_inserts:
    input:
        cbams = expand("combined_realign_all/{chrom}.Realigned.bam",chrom=config["chroms"]),
        cbams_bai = expand("combined_realign_all/{chrom}.Realigned.bam.bai",chrom=config["chroms"]),
        tsv = expand("combined_realign_all/{chrom}.Realigned.tsv",chrom=config["chroms"])
    output:
        bam="combined_realign_all/Realigned.tmp.bam",
        tsv="combined_realign_all/Realigned.tsv"
    threads: 1
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        depth=max_depth,
        ref=get_ref,
        bams=get_bam_list
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py merge --output-dir {params.base_dir}/combined_realign_all --tsv-prefix Realigned --bam-prefix Realigned --bam-list {params.bams}
        cd {params.base_dir}
        mv combined_realign_all/Realigned.bam combined_realign_all/Realigned.tmp.bam
        """

rule sort_realign_bams:
    input:
        bam="combined_realign_all/Realigned.tmp.bam"
    output:
        bam="combined_realign_all/Realigned.bam",
    threads: 1
    params:
        memory_per_thread="24G"
    shell:
        """
        samtools sort {input.bam} > {output.bam}
        rm {input.bam}
        """

rule classify_inserts:
    input:
        bam="combined_realign_all/Realigned.bam",
        bai="combined_realign_all/Realigned.bam.bai",
        realign_tsv="combined_realign_all/Realigned.tsv"
    output:
        tsv=temp("combined_realign_all/Realigned_classified.tsv")
    threads: 1
    params:
        memory_per_thread="64G",
        base_dir=get_base_dir,
        repbase=get_repbase,
        tsv_list=get_tsv_list,
        fastq_list=get_fastq_list,
        sample_list=get_sample_list
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py classify --sample-list {params.sample_list} --tsv-list {params.tsv_list} --realign-tsv {params.base_dir}/{input.realign_tsv} --annotation-file {params.repbase} --fastq-list {params.fastq_list} --output-tsv {params.base_dir}/{output.tsv}
        cd {params.base_dir}
        """

rule filter_inserts_chrom:
    input:
        tsv="combined_realign_all/Realigned_classified.tsv",
        bam="combined_realign_all/Realigned.bam",
        bai="combined_realign_all/Realigned.bam.bai"
    output:
        tsv="combined_realign_all/Realigned_classified_filtered.{chrom}.tsv"
    threads: 10
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        ref=get_ref,
        repbase=get_repbase,
        centromeres=get_centromeres,
        telomeres=get_telomeres,
        fastq_list=get_fastq_list
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py filter --threads {threads} --input-tsv {params.base_dir}/{input.tsv} --bam {params.base_dir}/{input.bam} --fastq-list {params.fastq_list} --reference-genome {params.ref} --output-tsv {params.base_dir}/{output.tsv} --centromeres {params.centromeres} --chromosome-list {wildcards.chrom}
        cd {params.base_dir}
        """

rule merge_filter_inserts:
    input:
        expand("combined_realign_all/Realigned_classified_filtered.{chrom}.tsv",chrom=config["chroms"])
    output:
        "combined_realign_all/Realigned_classified_filtered.tsv"
    threads: 1
    params:
        memory_per_thread="10G"
    shell:
        """
        head -1 combined_realign_all/Realigned_classified_filtered.chr1.tsv > {output}
        for file in {input}
        do
          tail -n+2 "$file" >> {output}
        done
        """

#rule summarize_data:
#    input:
#        realign_tsv="{sample}/Realigned_classified_filtered_{sample}.tsv",
#        tldr_table="{sample}/tldr_{sample}./{sample}.table.txt",
#    output:
#        tsv="simulation_results_{sample}.txt",
#        tsv_500="simulation_results_500_{sample}.txt",
#        tsv_2000="simulation_results_2000_{sample}.txt",
#        tsv_6000="simulation_results_6000_{sample}.txt",
#        tsv_LINE="simulation_results_LINE_{sample}.txt",
#        tsv_Alu="simulation_results_Alu_{sample}.txt",
#        tsv_SVA="simulation_results_SVA_{sample}.txt"
#    threads: 1
#    params:
#        memory_per_thread="64G",
#        somrit_test_script = srcdir("../scripts/get_simulation_metrics.py")
#    shell:
#        """
        #python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {input.xtea_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} > {output.tsv}
        #python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {input.xtea_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --max-size 500 > {output.tsv_500}
        #python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {input.xtea_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --min-size 500 --max-size 2000 > {output.tsv_2000}
        #python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {input.xtea_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --min-size 2000 > {output.tsv_6000}
        #python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {input.xtea_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --family LINE > {output.tsv_LINE}
        #python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {input.xtea_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --family Alu > {output.tsv_Alu}
        #python {params.somrit_test_script} --inserts-tsv {input.inserts_tsv} --spiked-reads {input.spike_in_reads} --tldr {input.tldr_table} --xtea {input.xtea_table} --somrit {input.realign_tsv} --sniffles {input.sniffles_vcf} --rep {wildcards.rep} --family SVA > {output.tsv_SVA}
#        """