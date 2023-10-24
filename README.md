# somrit-test

This repository contains scripts and snakemake rules to run the evaluation and analysis presented in the somrit pre-print [here](https://www.biorxiv.org/content/10.1101/2023.08.06.552193v2). 
The anlaysis is broken up into two sections. The first requires data from the [HPRC](https://humanpangenome.org/) and the second uses publically available data generated by [Gerdes et al 2022](https://www.nature.com/articles/s41467-022-35180-x)

The following tools are required. For some tools the path to the install folder must be specified in the ```project_config.yaml``` file. See the config file and/or the sections below to see which tools are required for which analysis. 

- [somrit](https://github.com/adcosta17/somrit)
- [tldr](https://github.com/adamewing/tldr)
- [xTea-Long](https://github.com/parklab/xTea)
- [sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [pbsim2](https://github.com/yukiteruono/pbsim2)
- python and matplotlib

## HPRC based analysis 

We used publically available data from the Human Pangenome Reference Consortium <> for the majority of our anlaysis. The HPRC has sequenced multiple karyotypically normal samples using Oxford Nanopore, Illumina and PacBio among other sequencing technologies. For HPRC samples HG00438, HG00621, HG00673, HG00735 and HG00741 we downloaded the following, replacing the SAMPLE placeholder with the specific sample ID: 

- Hifiasm based Diploid Assmemblies from: ```https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/SAMPLE/assemblies/year1_f1_assembly_v2_genbank/```
- Oxford Nanopore read sets from : ```https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/SAMPLE/raw_data/nanopore/```
- The RepeatMaster annotation of the diplid assembly from: ```https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/SAMPLE/assemblies/year1_f1_assembly_v2_genbank/annotation/repeat_masker/```

After downloading the ONT read sets were merged together into a single fastq while the maternal and paternal contig sets were renamed to mat.fa and pat.fa repsectively. All the HPRC data was downloaded into one top level folder, with a sub-folder per sample. Each sample's sub-folder had additional folders where fastqs were stored (fastq) and assembly contigs and annotations were stored (assembly). 

In addition to the data listed above for each sample we downloaded the RepeatMaster annotation of the reference genome GRCh38 from [here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/resources/hg38.rmsk.txt.gz)

### Analysis of Polymorphic Variation and Translocation Calls

To analyse somrit, tldr and xTea-Long's ability to detect polymorphic insertions we ran each sample individually. We first computed the total coverage of the ONT read set for the sample, using the coverage to identify the fraction of the read set required to downample to set fastq levels between 3x and 33x, listing these as the 'fractions' in the config file for the sample. We use the following 'fractions' for each HPRC sample. :

- HG00438: ```0.067,0.133,0.2,0.266,0.333,0.4,0.466,0.533,0.6,0.666,0.733```
- HG00621: ```0.0857,0.1714,0.2571,0.343,0.4286,0.5143,0.6,0.6857,0.7714,0.8571,0.9429```
- HG00673: ```0.0789,0.1579,0.2368,0.3158,0.3947,0.4736,0.5526,0.6316,0.7105,0.7895,0.8684```
- HG00735: ```0.088,0.176,0.265,0.353,0.441,0.529,0.618,0.706,0.794,0.882,0.971```
- HG00741: ```0.063,0.125,0.188,0.25,0.313,0.375,0.474,0.5,0.563,0.625,0.688```

We then ran the following snakemake rules, repeating it for each of the 5 samples. Each snakemake run for a sample should be run within a sample specific sub-folder. 

```
# First Run to generate the data
snakemake -s test.smk --configfile project_config.yaml all_sniffles_cutesv_somrit_xtea_tldr
# Second to summarize the data
snakemake -s test.smk --configfile project_config.yaml all_metrics
```

Lastly once the following has been run for all 5 samples we can run the following script to generate the plots shown in the paper. 
```
python scripts/plots.py --input-dir INPUT_DIR --somrit-dir SOMRIT_DIR --output-prefix OUTPUT_PREFIX --output-dir OUTPUT_DIR --samples SAMPLES --fractions FRACTIONS --coverage COVERAGE
```
INPUT_DIR and SOMRIT_DIR should be the same and should point to the top level folder that contains the sub-folders for each snakemake run. SAMPLES lists the sample names as a csv list. COVERAGE the coverage per sample read set in the same order as SAMPLES as a csv list. FRACTIONS is a list of lists. Each sample has a list of fraction values computed and set above to sepcifiy the downsampled coverage levels. These are represented as a csv list per sample, with each samples list seperated by a ':'. The call used to generate the plots shown in the paper from the 5 HPRC sample listed is shown below. 

```
python scripts/plots.py --input-dir . --somrit-dir . --output-prefix combined --output-dir plots/ --samples HG00438,HG00621,HG00673,HG00735,HG00741 --fractions 0.067,0.133,0.2,0.266,0.333,0.4,0.466,0.533,0.6,0.666,0.733:0.0857,0.1714,0.2571,0.343,0.4286,0.5143,0.6,0.6857,0.7714,0.8571,0.9429:0.0789,0.1579,0.2368,0.3158,0.3947,0.4736,0.5526,0.6316,0.7105,0.7895,0.8684:0.088,0.176,0.265,0.353,0.441,0.529,0.618,0.706,0.794,0.882,0.971:0.063,0.125,0.188,0.250,0.313,0.375,0.474,0.500,0.563,0.625,0.688 --coverage 45,35,38,34,48
```

### Analysis of Simulated Somatic Insertions

We used the same HPRC samples to generate simulated somatic insertions and evaulate somrit, tldr, xTea-Long and Sniffles2's performance of detection for simulated somatic insertions. 

We modified the number of replicates used from 3 to 12. In addition to the parameters specified above the simulated analysis requires the additional information be provided in the config files per sample

- ```mat_fa``` and ```pat_fa```: The paths to the maternal and paternal assembly contig sets respectively
- ```mat_rm``` and ```pat_rm```: The paths to the maternal and paternal assembly RepeatMasker annotations
- ```assembly_folder```: The top level folder where the HPRC data is stored
- ```pbsim_model```: The path to the pbsim2 module used for the pbsim2 simulation of the reads. Used R94.model for this analysis
- ```chrom_lengths```: The lengths of the GRCh38 chromosomes. Used utils/GRCh38_chromosome_sizes.tsv

We ran the following snakemake rules:

```
# To generate the data
snakemake -s simulation_test.smk --configfile project_config.yaml all_realign_tldr
# To summarize the data
snakemake -s simulation_test.smk --configfile project_config.yaml all
```

Once data was summarized we copied the summarized result files, prefixed with "simulation_results_" to sub-folder per sample before generating plots using the script listed below:

```
python scripts/simulation_plots.py --input-dir INPUT_DIR --output-prefix OUTPUT_PREFIX --output-dir OUTPUT_DIR --samples SAMPLES
```
INPUT_DIR should be the top level folder where each sample specific sub-folder containing the summarized results files is. SAMPLES is a csv list of samples used. 

## Gerdes et al based Analysis

In Addition to HPRC based evaluations we used publically available data from a recently published paper by [Gerdes et al 2022](https://www.nature.com/articles/s41467-022-35180-x) . In this paper the authors generated 5 HeLa cell line samples that contained novel insertions of a modified L1-mCherry construct. We downloaded used this data to compare somrit and tldr's ablility to identify novel insertion events. Each read set was placed in a sample specificic subfolder of the form "sample/fastq/sample.fastq.gz".

### L1-mCherry Insertion detection in treated HeLa Cells

We used the ```hela_project_config.yaml``` file and ran the following snakemake rules: 

```
snakemake -s hela_test.smk --configfile hela_project_config.yaml all
```

We then ran the following script to generate the comparison results

```
python scripts/get_hela_metrics.py --tldr TLDR --somrit SOMRIT --line-seq LINE_SEQ --l1mcherry-seq L1MCHERRY_SEQ --mcherry-seq MCHERRY_SEQ --insertions INSERTIONS
```
We used the following in our analysis TLDR is the ```combined/combined.tldr.table.txt``` and SOMRIT the ```combined_realign_all/Realigned_classified_filtered.tsv```. LINE_SEQ is ```utils/mouse_l1.fa```, L1MCHERRY_SEQ is ```utils/L1_EF1alpha-mCherry_spliced.fa```, MCHERRY_SEQ is ```utils/mCherry.fa``` and INSERTIONS is the truth set ```utils/mouse_truth.csv```.
