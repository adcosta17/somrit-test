#
# Data & Fastq related rules
#

def get_ref(wildcards):
    return config["reference"]

def get_fastq(wildcards):
    return config["fastq"]

def get_sample(wildcards):
    return config["full_sample"]

def get_chrom_list(wildcards):
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def get_base_dir(wildcards):
    return config["base_dir"]

def get_maternal(wildcards):
    return config["mat_fa"]

def get_paternal(wildcards):
    return config["pat_fa"]

def get_maternal_rm(wildcards):
    return config["mat_rm"]

def get_paternal_rm(wildcards):
    return config["pat_rm"]


rule sample_fastqs:
    output:
        fastq="HG{sample}_{frac}_{rep}.fastq.bgz"
    threads: 1
    log:
        "logs/fastq/HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/fastq/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="24G",
        fastq=get_fastq
    shell:
        """
        seqtk sample -s {wildcards.rep} {params.fastq} {wildcards.frac} | bgzip -c > {output.fastq} 2> {log}
        """

rule map_fastqs:
    input:
        fastq="HG{sample}_{frac}_{rep}.fastq.bgz"
    output:
        bam="HG{sample}_{frac}_{rep}.bam"
    threads: 10
    log:
        "logs/minimap2/HG{sample}_{frac}_{rep}.log"
    benchmark:
        "benchmark/minimap2/HG{sample}_{frac}_{rep}.txt"
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.ref} {input.fastq} | samtools sort > {output.bam} 2> {log}
        """


rule align_maternal:
    output:
        "HG{sample}.mat.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_maternal
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax asm5  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """

rule align_paternal:
    output:
        "HG{sample}.pat.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use=get_ref,
        fasta=get_paternal
    threads: 20
    shell:
        """
        minimap2 -t {threads} -ax asm5  {params.ref_to_use} {params.fasta} | samtools sort -o {output}
        """


rule get_rm_truth:
    input:
        pat="HG{sample}.pat.sorted.bam",
        mat="HG{sample}.mat.sorted.bam",
        pat_bai="HG{sample}.pat.sorted.bam.bai",
        mat_bai="HG{sample}.mat.sorted.bam.bai"
    output:
        "HG{sample}.rm_annotated.assembly_inserts_list.txt"
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/generate_truth_data.py"),
        chr_list=get_chrom_list,
        sample=get_sample,
        input_path=get_base_dir,
        pat_fasta=get_paternal,
        mat_fasta=get_maternal,
        pat_rm=get_paternal_rm,
        mat_rm=get_maternal_rm
    threads: 1
    shell:
        """
        python {params.script} --mat-bam {input.mat} --pat-bam {input.pat} --mat-fasta {params.mat_fasta} --pat-fasta {params.pat_fasta} --repeat-master-mat {params.mat_rm} --repeat-master-pat {params.pat_rm} > {output}
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
        "HG{sample}.mat.subset.bam"
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
        "HG{sample}.pat.subset.bam"
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
        pat_bam_index="HG{sample}.pat.subset.bam.bai",
        mat_fa="HG{sample}.mat.subset.fa",
        pat_fa="HG{sample}.pat.subset.fa"
    output:
        mat="HG{sample}.mat.insertions.tsv",
        pat="HG{sample}.pat.insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="36G"
    shell:
        """
        python {params.candidate_insertion_script} --bam {input.mat_bam} --fasta {input.mat_fa} --min-insertion-length 50 --min-mapq 20 --min-detected-inclusion-length 50 > {output.mat}
        python {params.candidate_insertion_script} --bam {input.pat_bam} --fasta {input.pat_fa} --min-insertion-length 50 --min-mapq 20 --min-detected-inclusion-length 50 > {output.pat}
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

rule summarize_time:
    output:
        benchmark="Time_Mem_HG{sample}_{frac}_{rep}.txt"
    threads: 1
    params:
        memory_per_thread="64G",
        xtea_before_txt="HG{sample}_{frac}_{rep}/classified_results.txt.Combined.txt",
        benchmark_script = srcdir("../scripts/summarize_benchmark.py"),
        tsv="Realigned_classified_filtered_HG{sample}_{frac}_{rep}.tsv",
        sniffles_before_vcf="Sniffles_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.vcf",
        tldr_before_txt="tldr_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.table.txt",
        cute_sv_before_vcf="CuteSV_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.vcf"
    shell:
        """
        python {params.benchmark_script} --sample {wildcards.sample} --fraction {wildcards.frac} --rep {wildcards.rep} > {output.benchmark}
        """

rule summarize_somrit:
    input:
       "HG{sample}.rm_annotated.assembly_inserts_list.txt"
    output:
        somrit_before_after="Somrit_xtea_tldr_before_after_HG{sample}_{frac}_{rep}.txt",
        somrit_before_after_all="Somrit_xtea_tldr_before_after_all_HG{sample}_{frac}_{rep}.txt"
    threads: 1
    params:
        memory_per_thread="64G",
        xtea_before_txt="HG{sample}_{frac}_{rep}/classified_results.txt",
        somrit_test_script = srcdir("../scripts/get_xtea_tldr_somrit_metrics.py"),
        tsv="Realigned_classified_filtered_HG{sample}_{frac}_{rep}.tsv",
        tldr_before_txt="HG{sample}_{frac}_{rep}.table.txt",
    shell:
        """
        python {params.somrit_test_script} --truth {input} --tldr-before {params.tldr_before_txt} --xtea-before {params.xtea_before_txt} --somrit {params.tsv} --fraction {wildcards.frac} --rep {wildcards.rep} --out-all {output.somrit_before_after_all} --out-rt {output.somrit_before_after}
        """

rule summarize_sniffles_cute_sv:
    input:
        mat_indels="HG{sample}.mat.insertions.tsv",
        pat_indels="HG{sample}.pat.insertions.tsv",
        mat_repbase="HG{sample}.mat.insertions.repbase_annotated.tsv",
        pat_repbase="HG{sample}.pat.insertions.repbase_annotated.tsv"
    output:
        fp_translocations="HG{sample}_Summary_translocations_{frac}_{rep}.txt",
        insertions_before_after="Sniffles_CuteSV_before_after_HG{sample}_{frac}_{rep}.txt",
        insertions_rt_before_after="RTs_Sniffles_CuteSV_before_after_HG{sample}_{frac}_{rep}.txt",
        diff_before_after="Diff_Sniffles_CuteSV_before_after_HG{sample}_{frac}_{rep}.txt",
        rt_diff_before_after="RT_Diff_Sniffles_CuteSV_before_after_HG{sample}_{frac}_{rep}.txt",
    threads: 1
    params:
        memory_per_thread="64G",
        sniffles_script = srcdir("../scripts/get_sniffles_cute_sv_inserts.py"),
        diff_script = srcdir("../scripts/get_sniffles_cute_sv_diff.py"),
        base_bam="HG{sample}_{frac}_{rep}.bam",
        base_bai="HG{sample}_{frac}_{rep}.bam.bai",
        realigned_bam="Realigned_HG{sample}_{frac}_{rep}.bam",
        realigned_bai="Realigned_HG{sample}_{frac}_{rep}.bam.bai",
        sniffles_before_vcf="Sniffles_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.vcf",
        sniffles_after_vcf="Sniffles_Realigned_HG{sample}_{frac}_{rep}/Realigned_HG{sample}_{frac}_{rep}.vcf",
        cute_sv_before_vcf="CuteSV_HG{sample}_{frac}_{rep}/HG{sample}_{frac}_{rep}.vcf",
        cute_sv_after_vcf="CuteSV_Realigned_HG{sample}_{frac}_{rep}/Realigned_HG{sample}_{frac}_{rep}.vcf",
    shell:
        """
        python {params.sniffles_script} --truth-mat {input.mat_indels} --truth-pat {input.pat_indels} --sniffles-before {params.sniffles_before_vcf} --sniffles-after {params.sniffles_after_vcf} --cutesv-before {params.cute_sv_before_vcf} --cutesv-after {params.cute_sv_after_vcf} --fraction {wildcards.frac} --rep {wildcards.rep} --insert-stats {output.insertions_before_after} --rt-insert-stats {output.insertions_rt_before_after} --translocation-stats {output.fp_translocations}
        python {params.diff_script} --truth-mat {input.mat_indels} --truth-pat {input.pat_indels} --sniffles-before {params.sniffles_before_vcf} --sniffles-after {params.sniffles_after_vcf} --cutesv-before {params.cute_sv_before_vcf} --cutesv-after {params.cute_sv_after_vcf} --fraction {wildcards.frac} --rep {wildcards.rep} --insert-stats {output.diff_before_after} --rt-insert-stats {output.rt_diff_before_after} --base-bam {params.base_bam} --realigned-bam {params.realigned_bam}
        """
