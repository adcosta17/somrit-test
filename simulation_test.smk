import os
import glob

rule all:
    input:
        expand("simulation_results_HG{s}.spike_in.{r}.txt",s=config["samples"],r=config["replicates"])

rule all_insert_seqs:
    input:
        expand("HG{s}/HG{s}.insert_added_sequences.tsv",s=config["samples"])

rule all_pbsim:
    input:
        expand("HG{s}/HG{s}.simulated_fastq_list.tsv",s=config["samples"])

rule all_insert_reads_fastq:
    input:
        expand("HG{s}/HG{s}.spanning_reads_list.tsv",s=config["samples"])

rule all_spike_in_fastq:
    input:
        expand("HG{s}/HG{s}.spike_in.{r}.fastq.gz",s=config["samples"],r=config["replicates"])

rule all_realign:
    input:
        expand("HG{s}/Realigned_HG{s}.spike_in.{r}.bam",s=config["samples"],r=config["replicates"])

rule all_realign_filtered:
    input:
        expand("HG{s}/Realigned_classified_filtered_HG{s}.spike_in.{r}.tsv",s=config["samples"],r=config["replicates"])

rule all_realign_tldr:
    input:
        expand("HG{s}/Realigned_classified_filtered_HG{s}.spike_in.{r}.tsv",s=config["samples"],r=config["replicates"]),
        expand("HG{s}/tldr_HG{s}.spike_in.{r}/HG{s}.spike_in.{r}.table.txt",s=config["samples"],r=config["replicates"]),
        expand("HG{s}/Sniffles_HG{s}.spike_in.{r}/HG{s}.spike_in.{r}.vcf",s=config["samples"],r=config["replicates"]),
        expand("HG{s}/HG{s}.spike_in.{r}/HG{s}.spike_in.{r}/classified_results.txt.Combined.txt",s=config["samples"],r=config["replicates"])

#include: "rules/data.smk"
#include: "rules/somrit.smk"
#include: "rules/mapping.smk"
#include: "rules/xtea.smk"
#include: "rules/tldr.smk"
#include: "rules/sniffles.smk"
#include: "rules/cuteSV.smk"
include: "rules/simulation.smk"
#include: "rules/deletion_validation.smk"
