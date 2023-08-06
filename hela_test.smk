import os
import glob

rule all:
    input:
        "combined_realign_all/Realigned_classified_filtered.tsv",
        "combined/combined.tldr.table.txt"



include: "rules/hela.smk"

