# Polarization

## Overview

This repository contains a python script to polorize a vcf which includes a target population and two outgroup populations.

## How to take the polarization decision ?

<img width="1890" height="1005" alt="image" src="https://github.com/user-attachments/assets/886c09e8-f9d1-42cc-8ec0-32247e44c4ae" />

## Dependencies

- python==3.12
- pysam==0.23.3
- tqdm==4.67.1

## Running the script

Here is an exemple of how to use this script.

### Command line

python Polarisation.py -i <VCF> -t <OUT_TABLE> --pop <pop.txt> --out1 <out1.txt> --out2 <out2.txt>

### Arguments

-i <VCF> Name of the VCF containing the target population (pop.txt) and the two outgroup populations (out1.txt and out2.txt). 
-t <OUT_TABLE> Name of the output table countaing the following columns :
  "CHROM   POS     REF_ORIG        ALT_ORIG        ANC_ALLELE      DER_ALLELE      CASE    STATUS"
  '_ORIG' refers to the information at this locus in the vcf before the polarization. 'ANC_ALLELE' and 'DER_ALLELE' refer respectively to the the ancestral site and the      derived allele after the polarization. 
  'CASE' and 'STATUS' inform us on wheter the locus could be polarized. If it couldn't (due to missing data or insolvable alleles distribution), 'STATUS' will display the    string 'NOPOL'. Otherwise, the string will inform on which case scenario of allele distribution allowed us to polarize, along with the column 'CASE'. Note that the case    scenario rarely matters to the user, in which case it's only of matter of wheter the column 'STATUS' displays "NOPOL" (no polarization at this locus) or not (polarization at this locus).
--pop Name of the file countaining the names of all the individuals in the target population included in the VCF, with one individual name per line.
--out1.txt Name of the file countaining the names of all the individuals in the first outgroup population included in the VCF, with one individual name per line.
--out2.txt Name of the file countaining the names of all the individuals in the second outgroup population included in the VCF, with one individual name per line.

  
