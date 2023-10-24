import os
import glob

rule all:
    input:
        "combined_realign_all/Realigned_classified_filtered.tsv",
        "combined/combined.tldr.table.txt"

rule all_realign:
    input:
        "combined_realign_all/Realigned_classified_filtered.tsv"

rule all_normalized:
    input:
        expand("{s}/{s}_normalized.txt", s=config["samples"])

include: "rules/mice.smk"

