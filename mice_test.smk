import os
import glob

rule all:
    input:
        "combined_realign_all/Realigned_classified_filtered.tsv",
        "combined/combined.tldr.table.txt"

rule all_realign:
    input:
        "combined_realign_all/Realigned_classified_filtered.tsv"


include: "rules/mice.smk"

