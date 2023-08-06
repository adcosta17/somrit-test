import os
import glob

rule all_pbsim:
    input:
        "simulated_fastq_list.tsv"


#include: "rules/data.smk"
#include: "rules/somrit.smk"
#include: "rules/mapping.smk"
#include: "rules/xtea.smk"
#include: "rules/tldr.smk"
#include: "rules/sniffles.smk"
#include: "rules/cuteSV.smk"
include: "rules/nested_inserts.smk"
#include: "rules/deletion_validation.smk"
