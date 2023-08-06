
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
    return int(float(wildcards.frac)*600)

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
        temp("HG{sample}.mat.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_maternal
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule align_paternal:
    output:
        temp("HG{sample}.pat.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_paternal
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax asm20  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """


rule get_subset_assemblies:
    input:
        pat="HG{sample}.pat.sorted.bam",
        mat="HG{sample}.mat.sorted.bam",
        pat_bai="HG{sample}.pat.sorted.bam.bai",
        mat_bai="HG{sample}.mat.sorted.bam.bai"
    output:
        mat="HG{sample}.mat.subset.fa",
        pat="HG{sample}.pat.subset.fa"
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_aligned_asssembly.py"),
        chr_list=get_chrom_list,
        sample=get_sample,
        input_path=get_base_dir
    threads: 1
    shell:
        """
        python {params.script} --chrom-list {params.chr_list} --full-sample {params.sample} --sample HG{wildcards.sample} --input-path {params.input_path}
        """


rule map_mat_subset_to_ref:
    input:
        mat="HG{sample}.mat.subset.fa"
    output:
        temp("HG{sample}.mat.subset.bam")
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {params.ref_to_use} {input.mat} | samtools sort -o {output}
        """

rule map_pat_subset_to_ref:
    input:
        pat="HG{sample}.pat.subset.fa"
    output:
        temp("HG{sample}.pat.subset.bam")
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref
    threads: 10
    shell:
        """
        minimap2 -t {threads} -ax asm20 {params.ref_to_use} {input.pat} | samtools sort -o {output}
        """

rule get_pat_and_mat_ref_insertions:
    input:
        mat_bam="HG{sample}.mat.subset.bam",
        mat_bam_index="HG{sample}.mat.subset.bam.bai",
        pat_bam="HG{sample}.pat.subset.bam",
        pat_bam_index="HG{sample}.pat.subset.bam.bai"
    output:
        mat="HG{sample}.mat.insertions.tsv",
        pat="HG{sample}.pat.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        """
        python {params.candidate_insertion_script} --bam {input.mat_bam} --min-insertion-length 50 --min-mapq 20 --min-detected-inclusion-length 50 > {output.mat}
        python {params.candidate_insertion_script} --bam {input.pat_bam} --min-insertion-length 50 --min-mapq 20 --min-detected-inclusion-length 50 > {output.pat}
        """

rule run_convert_pat_and_mat_ref_insertions_to_fasta:
    input:
        mat="HG{sample}.mat.insertions.tsv",
        pat="HG{sample}.pat.insertions.tsv"
    output:
        mat=temp("HG{sample}.mat.insertions.fa"),
        pat=temp("HG{sample}.pat.insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_assembly_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        """
        python {params.candidate_insertion_conversion_script} --input {input.mat} > {output.mat}
        python {params.candidate_insertion_conversion_script} --input {input.pat} > {output.pat}
        """

rule make_lastdb:
    input:
        config["repbase"]
    output:
        config["repbase"] + ".lastdb.suf"
    threads: 1
    params:
        memory_per_thread="16G"
    shell:
        "lastdb {config[repbase]}.lastdb {input}"


rule map_insertion_sequences_last:
    input:
        fa="{base}.fa",
        db=config["repbase"] + ".lastdb.suf"
    output:
        "{base}.mapped_to_repbase.last.tab"
    threads: 8
    params:
        memory_per_thread="12G"
    shell:
        "lastal -P {threads} -f tab -r1 -a1 -b1 {config[repbase]}.lastdb {input.fa} > {output}"

rule run_pat_and_mat_ref_assembly_insertion_annotation:
    input:
        mat="HG{sample}.mat.insertions.tsv",
        pat="HG{sample}.pat.insertions.tsv",
        mat_tab="HG{sample}.mat.insertions.mapped_to_repbase.last.tab",
        pat_tab="HG{sample}.pat.insertions.mapped_to_repbase.last.tab"
    output:
        mat="HG{sample}.mat.insertions.repbase_annotated.tsv",
        pat="HG{sample}.pat.insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_assembly_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        """
        python {params.candidate_insertion_annotation_script} --input {input.mat} --last {input.mat_tab} --min-mapped-fraction 0 > {output.mat}
        python {params.candidate_insertion_annotation_script} --input {input.pat} --last {input.pat_tab} --min-mapped-fraction 0 > {output.pat}
        """

rule zip_fastq:
    input:
        "{prefix}.fastq"
    output:
        "{prefix}.fastq.gz"
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "bgzip -i {input}"

rule map_fastqs_to_assembly:
    output:
        bam=temp("HG{sample}/HG{sample}.{hap}.bam")
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_assembly,
        fastq=get_fastq
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.ref} {params.fastq} | samtools sort > {output.bam}
        """


rule get_deletion_fastqs:
    input:
        mat_inserts="HG{sample}.mat.insertions.repbase_annotated.tsv",
        pat_inserts="HG{sample}.pat.insertions.repbase_annotated.tsv",
        reads_mat_bam="HG{sample}/HG{sample}.mat.bam",
        reads_pat_bam="HG{sample}/HG{sample}.pat.bam",
        reads_mat_bam_index="HG{sample}/HG{sample}.mat.bam.bai",
        reads_pat_bam_index="HG{sample}/HG{sample}.pat.bam.bai",
    output:
        tsv="HG{sample}/HG{sample}.ct.annotation.tsv",
        fc="HG{sample}/HG{sample}.fastq",
        reads="HG{sample}/HG{sample}.ct.reads_to_ignore.txt"
    threads: 1
    params:
        script=srcdir("../scripts/get_deletion_fastqs.py"),
        memory_per_thread="64G",
        fastq=get_fastq,
        sample=get_sample,
        centromeres=get_centromeres,
        telomeres=get_telomeres
    shell:
        """
        python {params.script} --output-fastq {output.fc} --reads-to-ignore {output.reads} --input-fastq {params.fastq} --sample {params.sample} --mat-ref-inserts {input.mat_inserts} --pat-ref-inserts {input.pat_inserts} --mat-bam {input.reads_mat_bam} --pat-bam {input.reads_pat_bam} --centromeres {params.centromeres} > {output.tsv}
        """


rule sample_deletion_fastqs:
    input:
        fastq="HG{sample}/HG{sample}.fastq.gz"
    output:
        fastq=temp("HG{sample}_deletion/HG{sample}_{frac}_{rep}.fastq.bgz")
    threads: 1
    params:
        memory_per_thread="24G",
        fastq="HG{sample}_deletion/HG{sample}_{frac}_{rep}.fastq",
        script=srcdir("../scripts/sample_fastqs.py")
    shell:
        """
        python {params.script} --input {input.fastq} --output {params.fastq} --fraction {wildcards.frac} --seed {wildcards.rep}
        bgzip -i {params.fastq}
        """

rule map_deletion_fastqs:
    input:
        fastq="HG{sample}_deletion/HG{sample}_{frac}_{rep}.fastq.bgz"
    output:
        bam=temp("HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam")
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.ref} {input.fastq} | samtools sort > {output.bam} 2> {log}
        """

rule xtea_deletion:
    input:
        bam="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam.bai"
    output:
        tsv="HG{sample}_deletion/HG{sample}_{frac}_{rep}/classified_results.txt.Combined.txt"
    threads: 10
    benchmark:
        "benchmark/xtea_deletion/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        echo HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} > HG{wildcards.sample}_deletion/HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_sample_id.txt
        echo HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} {input.bam} > HG{wildcards.sample}_deletion/HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_long_read_bam_list.txt
        {config[xtea_dir]}/bin/xtea_long -i HG{wildcards.sample}_deletion/HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_sample_id.txt -b HG{wildcards.sample}_deletion/HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_long_read_bam_list.txt -p . -o HG{wildcards.sample}_deletion/HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_submit_jobs.sh --clean --rmsk {config[xtea_dir]}/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r {params.ref} --cns {config[xtea_dir]}/consensus/LINE1.fa --rep {config[xtea_dir]} --xtea {config[xtea_dir]}/xtea_long/ -f 31 -y 15 -n 10 -m 190 --clean
        rm HG{wildcards.sample}_deletion/HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}_submit_jobs.sh
        cd HG{wildcards.sample}_deletion/HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}
        chmod +x run_xTEA_pipeline.sh
        cd ../
        ./HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}/run_xTEA_pipeline.sh
        """


rule sniffles_deletion:
    input:
        bam="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam.bai"
    output:
        vcf="HG{sample}_deletion/Sniffles_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.vcf"
    threads: 10
    log:
        "logs/sniffles_deletion/HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/sniffles_deletion/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --non-germline &> {log}
        """

rule tldr_deletion:
    input:
        bam="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam.bai"
    output:
        tsv="HG{sample}_deletion/tldr_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.table.txt"
    threads: 10
    benchmark:
        "benchmark/tldr_deletion/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        cd HG{wildcards.sample}_deletion/tldr_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep}
        {config[tldr_dir]}/tldr/tldr -b ../{input.bam} -r {params.ref} -e {config[tldr_dir]}/ref/teref.ont.human.fa -p {threads}
        """


rule extract_inserts_deletion_fastqs:
    input:
        bam="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam.bai",
        fastq="HG{sample}_deletion/HG{sample}_{frac}_{rep}.fastq.bgz"
    output:
        tsv="HG{sample}_deletion/HG{sample}_{frac}_{rep}.tsv",
        merged=temp("HG{sample}_deletion/HG{sample}_{frac}_{rep}.merged.txt")
    threads: 10
    log:
        "logs/somrit_deletion/HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit_deletion/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="8G",
        base_dir=get_base_dir
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-merged {params.base_dir}/{output.merged} --output-tsv {params.base_dir}/{output.tsv} --fastq-file {params.base_dir}/{input.fastq} --threads {threads} &> {params.base_dir}/{log}
        cd {params.base_dir}
        """

rule realign_inserts_deletion_fastqs:
    input:
        bam="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam",
        tsv="HG{sample}_deletion/HG{sample}_{frac}_{rep}.tsv",
        fastq="HG{sample}_deletion/HG{sample}_{frac}_{rep}.fastq.bgz"
    output:
        bam="HG{sample}_deletion/{chrom}.Realigned_HG{sample}_{frac}_{rep}.bam",
        tsv="HG{sample}_deletion/{chrom}.Realigned_HG{sample}_{frac}_{rep}.tsv"
    threads: 10
    log:
        "logs/somrit_deletion/{chrom}.Realigned_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit_deletion/{chrom}.Realigned_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        depth=max_depth,
        ref=get_ref
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.base_dir}/{input.bam} --tsv-list {params.base_dir}/{input.tsv} --fastq-list {params.base_dir}/{input.fastq} --output-dir {params.base_dir}/HG{wildcards.sample}_deletion --tsv-prefix Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --bam-prefix Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 10000 --max-depth {params.depth} --chromosome-list {wildcards.chrom} --only-realign &> {params.base_dir}/{log}
        cd {params.base_dir}
        """


rule merge_realign_inserts_deletion_fastqs:
    input:
        bam="HG{sample}_deletion/HG{sample}_{frac}_{rep}.bam",
        cbams = expand("HG{{sample}}_deletion/{chrom}.Realigned_HG{{sample}}_{{frac}}_{{rep}}.bam",chrom=config["chroms"]),
        tsv = expand("HG{{sample}}_deletion/{chrom}.Realigned_HG{{sample}}_{{frac}}_{{rep}}.tsv",chrom=config["chroms"])
    output:
        bam=temp("HG{sample}_deletion/Realigned_HG{sample}_{frac}_{rep}.bam"),
        tsv="HG{sample}_deletion/Realigned_HG{sample}_{frac}_{rep}.tsv"
    threads: 1
    log:
        "logs/somrit_deletion/Realigned_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit_deletion/Realigned_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        depth=max_depth,
        ref=get_ref
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py merge --output-dir {params.base_dir}/HG{wildcards.sample}_deletion --bam-list {params.base_dir}/{input.bam} --tsv-prefix Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --bam-prefix Realigned_HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} &> {params.base_dir}/{log}
        cd {params.base_dir}
        """

rule classify_inserts_deletion_fastqs:
    input:
        tsv="HG{sample}_deletion/HG{sample}_{frac}_{rep}.tsv",
        bam="HG{sample}_deletion/Realigned_HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_deletion/Realigned_HG{sample}_{frac}_{rep}.bam.bai",
        realign_tsv="HG{sample}_deletion/Realigned_HG{sample}_{frac}_{rep}.tsv",
        fastq="HG{sample}_deletion/HG{sample}_{frac}_{rep}.fastq.bgz"
    output:
        tsv=temp("HG{sample}_deletion/Realigned_classified_HG{sample}_{frac}_{rep}.tsv")
    threads: 1
    log:
        "logs/somrit_deletion/Realigned_classified_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit_deletion/Realigned_classified_HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        repbase=get_repbase
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py classify --bam-list {params.base_dir}/{input.bam} --sample-list HG{wildcards.sample}_{wildcards.frac}_{wildcards.rep} --tsv-list {params.base_dir}/{input.tsv} --realign-tsv {params.base_dir}/{input.realign_tsv} --annotation-file {params.repbase} --fastq-list {params.base_dir}/{input.fastq} --output-tsv {params.base_dir}/{output.tsv} &> {params.base_dir}/{log}
        cd {params.base_dir}
        """

rule filter_inserts_deletion_fastqs:
    input:
        tsv="HG{sample}_deletion/Realigned_classified_HG{sample}_{frac}_{rep}.tsv",
        bam="HG{sample}_deletion/Realigned_HG{sample}_{frac}_{rep}.bam",
        bai="HG{sample}_deletion/Realigned_HG{sample}_{frac}_{rep}.bam.bai",
        fastq="HG{sample}_deletion/HG{sample}_{frac}_{rep}.fastq.bgz"
    output:
        tsv="HG{sample}_deletion/Realigned_classified_filtered_HG{sample}_{frac}_{rep}.tsv"
    threads: 10
    log:
        "logs/somrit_deletion/Realigned_classified_filtered_HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/somrit_deletion/Realigned_classified_filtered_HG{sample}_{frac}_{rep}.txt"
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
        python somrit.py filter --threads {threads} --input-tsv {params.base_dir}/{input.tsv} --bam {params.base_dir}/{input.bam} --fastq-list {params.base_dir}/{input.fastq} --reference-genome {params.ref} --centromeres {params.centromeres} --telomeres {params.telomeres} --output-tsv {params.base_dir}/{output.tsv} &> {params.base_dir}/{log}
        cd {params.base_dir}
        """



rule summarize_deletion_data:
    input:
        tsv="HG{sample}_deletion/Realigned_classified_filtered_HG{sample}_{frac}_{rep}.tsv",
        sniffles_vcf="HG{sample}_deletion/Sniffles_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.vcf",
        tldr_txt="HG{sample}_deletion/tldr_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.table.txt",
        xtea_txt="HG{sample}_deletion/HG{sample}_{frac}_{rep}/classified_results.txt.Combined.txt",
        mat_repbase="HG{sample}.mat.insertions.repbase_annotated.tsv",
        pat_repbase="HG{sample}.pat.insertions.repbase_annotated.tsv"
    output:
        somrit_before_after="Somrit_xtea_tldr_before_after_HG{sample}_{frac}_{rep}.txt"
    threads: 1
    params:
        memory_per_thread="64G",
        somrit_test_script = srcdir("../scripts/get_xtea_tldr_somrit_deletion_metrics.py")
    shell:
        """
        python {params.somrit_test_script} --truth-mat {input.mat_repbase} --truth-pat {input.pat_repbase} --tldr {input.tldr_txt} --xtea {input.xtea_txt} --somrit {input.tsv} --sniffles {input.sniffles_vcf} --fraction {wildcards.frac} --rep {wildcards.rep} > {output.somrit_before_after}
        """