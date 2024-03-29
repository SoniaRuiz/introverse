######################################################
## 1. FILE STRUCTURE
#######################################################

1. init.R
2. junction_pairing.R
3. database_SQL_helper.R
4. database_SQL_generation.R
 
 
######################################################
## 2. DEPENDENCIES
######################################################

/*****************************************************
** These are files stored within the
** './dependencies' folder
******************************************************/

1. 'CNC_CDTS_CONS_gr.rds'  
=======================
Contains the inter-species conservation scores across primates and the contraint regions across humans of every 10b-long region of the human genome.
Contact 'mina.ryten@ucl.ac.uk' for data download queries.


2. 'hg38-blacklist.v2.bed'
=======================
File that contains the ENCODE blacklist regions found in the hg38.
This is the paper: https://www.nature.com/articles/s41598-019-45839-z
This file has been downloaded from: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz 


3. 'major_introns_tidy.rds' and 'minor_introns_tidy.rds'
========================================================
Original BED files downloaded from the IAOD database:
https://introndb.lerner.ccf.org/sitepages/downloads/


4. 'MANE.GRCh38.v1.0.ensembl_genomic.gtf'
=========================================
MANE info: https://www.ncbi.nlm.nih.gov/refseq/MANE/
gtf file downloaded from: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz


5. 'MT_genes.rds'
==================
"ENSG00000210049" "ENSG00000211459" "ENSG00000210077" "ENSG00000210082" 
"ENSG00000209082" "ENSG00000198888" "ENSG00000210100"
"ENSG00000210107" "ENSG00000210112" "ENSG00000198763" "ENSG00000210117" 
"ENSG00000210127" "ENSG00000210135" "ENSG00000210140"
"ENSG00000210144" "ENSG00000198804" "ENSG00000210151" "ENSG00000210154" 
"ENSG00000198712" "ENSG00000210156" "ENSG00000228253"
"ENSG00000198899" "ENSG00000198938" "ENSG00000210164" "ENSG00000198840" 
"ENSG00000210174" "ENSG00000212907" "ENSG00000198886"
"ENSG00000210176" "ENSG00000210184" "ENSG00000210191" "ENSG00000198786" 
"ENSG00000198695" "ENSG00000210194" "ENSG00000198727"
"ENSG00000210195" "ENSG00000210196"


6. 'clinvar_intronic_tidy.rda'
==============================
.tsv data downloaded from https://www.ncbi.nlm.nih.gov/variation/view


7. 'Homo_sapiens.GRCh38.105.chr.gtf'
===================================
Data downloaded from: http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/


8. 'Homo_sapiens.GRCh38.dna.primary_assembly.fa', 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
===================================================================================================
## FASTA file
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --no-check-certificate
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx ./Homo_sapiens.GRCh38.dna.primary_assembly.fa


9. MAXENTSCAN score software
==============================
Software downloaded from 'http://hollywood.mit.edu/burgelab/software.html'
"../fordownload/..."

10. Bedtools
=============
bedtools software


######################################################
## 3. TESTING FRAMEWORK
######################################################

* Unit test implemented by: Guillermo Rocamora - https://github.com/guillermo1996
* Library used: https://testthat.r-lib.org/
