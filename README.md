# BENG 202 Project 
---
## Problem Statement: Identifying Potential Fluorophore-Binding Aptamers

The objective is to identify all potential mRNA k-mers that can bind to a specific fluorophore by analyzing both nucleotide similarity and structural similarity to known fluorophore-binding aptamers.

## Workplan

Many genomes encode structured RNAs that may have fluorophore-binding properties. To systematically identify such RNAs, we will integrate bioinformatics approaches with genomic and transcriptomic data analysis. Our workflow consists of four primary steps followed by computational analysis, ensuring a robust pipeline for discovering novel fluorophore-binding RNA motifs.

### Step 1 - Selecting a candidate organism

Before searching for fluorophore-binding sequences, we first determine an appropriate model organism. We choose Escherichia coli as our candidate because:

1. It is widely used in synthetic biology for engineering RNA-based tools.
2. Its transcriptome is well-annotated in public databases (e.g., Ensembl).
3. It is experimentally tractable for in vitro and in vivo validation.

### Step 2 - Retrieving transcriptomics data

We will look at genomic/transcriptomic datasets present in Ensembl BioMart to locate and download the latest transcriptomic mRNA data and download the transcript sequences in the FASTA format  

### Step 3 - Filtering candidate RNA sequences: We then narrow down our pool of mRNAs using the following steps:

Select  sequences that are annotated as non-coding RNAs
Extract previously reported riboswitches and aptamer-like sequences from this dataset

### Step 4 - Computational analysis of fluorophore binding sequences
Once we have obtained the transcriptomic dataset for Escherichia coli, we apply computational methods to identify k-mers that have the potential to bind fluorophores. The analysis involves evaluating each candidate RNA k-mer based on its sequence similarity and structural similarity to known fluorophore-binding aptamers.

#### Generating k-mers from RNA Sequences
Given k is the length of the known fluorophore binding sequences, we extract all possible k-mers (subsequences of length k) from the RNA transcripts.

#### Evaluation of nucleotide similarity (SeqSim function)
Each k-mer is compared to the set of known fluorophore-binding RNAs (F) using SeqSim(str,F) to assess the sequence similarity

#### Evaluation of structural similarity (StructSim function)
We predict the secondary structure of each k-mer using RNA folding algorithms
We compare the secondary structure representation (ssr) of each k-mer to the ssr of the known fluorophore-binding aptamers in F

#### Filtering candidate k-mers (Scoring function)
Each k-mer receives a final score based on a weighted combination of nucleotide similarity and structural similarity.
We apply a score threshold (t) to filter out low-scoring k-mers, retaining only those that are most likely to bind fluorophores.

#### Obtaining final k-mers
The final output is a set of high-scoring k-mers (H) that are strong candidates for fluorophore binding.
These k-mers can be further validated experimentally or used to design novel fluorophore-binding RNAs.
