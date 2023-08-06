import os
import glob

rule all_mat_pat_bam:
    input:
        expand("HG{s}.pat.sorted.bam.bai", s=config["samples"]),
        expand("HG{s}.mat.sorted.bam.bai", s=config["samples"])

rule all_sorted_bams:
    input:
        expand("HG{s}_{f}_{r}.bam", f=config["fractions"], r=config["replicates"], s=config["samples"])

rule all_realigned_bams:
    input:
        expand("Realigned_HG{s}_{f}_{r}.bam", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Realigned_HG{s}_{f}_{r}.bam.bai", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_xtea:
    input:
        expand("HG{s}_{f}_{r}/classified_results.txt.Combined.txt", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_tldr:
    input:
        expand("tldr_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.table.txt", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("tldr_HG{s}_{f}_{r}/HG{s}_{f}_{r}.table.txt", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_sniffles_cutesv:
    input:
        expand("Sniffles_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Sniffles_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_sniffles_cutesv_realigned:
    input:
        expand("Sniffles_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_somrit:
    input:
        expand("Realigned_classified_filtered_HG{s}_{f}_{r}.tsv", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_sniffles_cutesv_somrit:
    input:
        expand("Sniffles_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Sniffles_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Realigned_classified_filtered_HG{s}_{f}_{r}.tsv", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_sniffles_cutesv_somrit_xtea_tldr:
    input:
        expand("Sniffles_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Sniffles_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Realigned_classified_filtered_HG{s}_{f}_{r}.tsv", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("HG{s}_{f}_{r}/classified_results.txt.Combined.txt", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("tldr_HG{s}_{f}_{r}/HG{s}_{f}_{r}.table.txt", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_somrit_xtea_tldr:
    input:
        expand("Realigned_classified_filtered_HG{s}_{f}_{r}.tsv", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("HG{s}_{f}_{r}/classified_results.txt.Combined.txt", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("tldr_HG{s}_{f}_{r}/HG{s}_{f}_{r}.table.txt", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("HG{s}.rm_annotated.assembly_inserts_list.txt",s=config["samples"])

rule all_xtea_before:
    input:
        expand("HG{s}_{f}_{r}/classified_results.txt.Combined.txt", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_data:
    input:
        expand("Realigned_classified_filtered_HG{s}_{f}_{r}.tsv", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Sniffles_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Sniffles_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("tldr_HG{s}_{f}_{r}/HG{s}_{f}_{r}.table.txt", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("tldr_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.table.txt", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("HG{s}_{f}_{r}/classified_results.txt.Combined.txt", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Realigned_HG{s}_{f}_{r}/classified_results.txt.Combined.txt", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("HG{s}.mat.insertions.tsv",s=config["samples"]),
        expand("HG{s}.pat.insertions.tsv",s=config["samples"]),
        expand("HG{s}.mat.insertions.repbase_annotated.tsv",s=config["samples"]),
        expand("HG{s}.pat.insertions.repbase_annotated.tsv",s=config["samples"])

rule all_metrics:
    input:
        expand("HG{s}_Summary_translocations_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Sniffles_CuteSV_before_after_HG{s}_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("RTs_Sniffles_CuteSV_before_after_HG{s}_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Somrit_xtea_tldr_before_after_HG{s}_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Time_Mem_HG{s}_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_sniffles_cute_sv_metrics:
    input:
        expand("HG{s}_Summary_translocations_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Sniffles_CuteSV_before_after_HG{s}_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("RTs_Sniffles_CuteSV_before_after_HG{s}_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_sniffles_cutesv_somrit_tldr:
    input:
        expand("Sniffles_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Sniffles_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_HG{s}_{f}_{r}/HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("CuteSV_Realigned_HG{s}_{f}_{r}/Realigned_HG{s}_{f}_{r}.vcf", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Realigned_classified_filtered_HG{s}_{f}_{r}.tsv", f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("tldr_HG{s}_{f}_{r}/HG{s}_{f}_{r}.table.txt", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_fastq:
    input:
        expand("HG{s}_{f}_{r}.fastq.bgz", f=config["fractions"], r=config["replicates"],s=config["samples"])

rule all_somrit_metrics:
    input:
        expand("Somrit_xtea_tldr_before_after_HG{s}_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"]),
        expand("Time_Mem_HG{s}_{f}_{r}.txt",f=config["fractions"], r=config["replicates"],s=config["samples"])

include: "rules/data.smk"
include: "rules/somrit.smk"
include: "rules/mapping.smk"
include: "rules/xtea.smk"
include: "rules/tldr.smk"
include: "rules/sniffles.smk"
include: "rules/cuteSV.smk"
#include: "rules/deletion_validation.smk"
