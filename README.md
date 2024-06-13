# coprofilling-pipeline
Introduction
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
This workflow provides a comprehensive set of analyses for genomic data, including CNV, SNP, mitochondrial content, and fragments analysis, using a combination of command-line tools and custom Python scripts.

Install
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
python3

Aneufinder

freebayes

os

csv

math

pysam

numpy

Metrics

argparse

samtools

pipeline

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Step 1: Alignment using BWA
Command: bwa mem reference.fasta reads.fastq > aligned.sam
This command aligns sequencing reads to the reference genome using the BWA mem algorithm, producing an SAM file.
Step 2: Processing Mapped Data
Scripts:
getAllBulkFiles_from_mappedData.py: Processes bulk sequencing data to generate sorted and deduplicated files.
getReadsidBarcodeDic_and_seperateSam.py: Separates SAM files based on barcode information.
getBam_for_seperatedSam.py: Converts separated SAM files into BAM format.
Step 3: Preparing Single-Cell BAM Files
Script: copy_SCbams_inOneDir.py
This script collects individual cell BAM files into a single directory for further analysis with Aneufinder.
Step 4: Copy Number and Single-Nucleotide Variant Analysis
Copy Number Variation:
Tool: Aneufinder for detecting copy number variations (CNVs).
Script: get_CNV_statiscsInfo.py: Summarizes CNV statistics for each sample. For bulk samples, get_bulkCNV_statiscsInfo.py is used.
Single-Nucleotide Variation:
Tool: Freebayes for single-nucleotide variation detection.
Scripts:
getSNPstatisticInfo.py: Generates SNP statistics for each barcode at chromosome precision.
getSNPstatisticInfoEachVariant_reduced.py: Focuses on counting variant sites without detailed variant information.
Chromosome-specific scripts for humans and mice (getSNPstatisticInfoEachBin_forHuman.py, getSNPstatisticInfoEachBinWithChrom_forHuman.py, etc.) provide SNP counts within specific window lengths.
Step 4: Mitochondrial Content Analysis
Script: get_UMIRemoveDup_MTRatioInfo.py
Calculates the mitochondrial content ratio after removing duplicate UMIs.
Step 5: Fragments Information
Script: getUMItypes.py
Statistics on fragment information based on UMI types.
Step 6: Sequencing Depth Information
Script: get_depthResult_and_statiscsInfo.py
Provides statistics on sequencing depth.
Step 7: Sequencing Coverage Information
Script: get_coverageResult_and_statiscsInfo.py
Provides statistics on sequencing coverage.
Step 11: Nucleus Ratio Information
Script: get_UMIRemoveDup_CellNucleusRatioInfo.py
Calculates the nucleus ratio in cells after duplicate UMI removal.
Step 12: Fragments Per Bin
Script: getPerBinFragments.py
Statistics on the number of fragments per specified genomic window.
Step 13: Chromosome-specific Fragments Information
Script: getPerChromFragments.py
Provides the number of fragments per chromosome.
Step 14: Chromosomal Window-specific Fragments Information
Script: getPerBinWithChromFragments_GEMB.py
Generates statistics on fragments within specified windows including chromosome information.
Step 15: Subclone-specific SNP Information
Scripts:
getSNPstatisticInfoEachBinWithChrom_withSubclone_forHuman.py
getSNPstatisticInfoEachBinWithChrom_withSubclone_forMouse.py
These scripts generate SNP statistics for different subclones within specified genomic windows.
