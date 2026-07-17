# Polarization

## Overview

This repository includes a python script to polarize a vcf which contains a target population and two outgroup populations.

## How to take the polarization decision ?

<img width="1890" height="1005" alt="image" src="https://github.com/user-attachments/assets/886c09e8-f9d1-42cc-8ec0-32247e44c4ae" />

## Dependencies

# Dependencies
- python==3.12
- Check [`requirements.txt`](./requirements.txt) for needed python packages

## Installation

```bash
pip install -r requirements.txt

## Running the script

### Warnings

- The VCF output might be corrupted for some analysis, as it is built up manually in the script. This issue will be fixed hopefully soon, in the meantime please use the output dataframe.
- Please make sure that the VCF contains only bi-allelic SNPs. Adding monomorphic sites would only slow down the script, and it does not handle positions with more than two distinct alleles or indels. You can for example run the following command : bcftools view -m2 -M2 -v snps input.vcf.gz -Oz -o filtered.vcf.gz


### Command line

python Polarisation.py -i <VCF> -o <OUT_VCF> -t <OUT_TABLE> --pop <pop.txt> --out1 <out1.txt> --out2 <out2.txt>

### Arguments

-i <VCF> Name of the VCF containing the target population (pop.txt) and the two outgroup populations (out1.txt and out2.txt).

-o <OUT_VCF> Name of the output VCF after polarization (to use with caution)

-t <OUT_TABLE> Name of the output table containg the following columns :

"CHROM   POS     REF_ORIG        ALT_ORIG        ANC_ALLELE      DER_ALLELE      CASE    STATUS"

'_ORIG' refers to the information at this locus in the vcf before the polarization. 'ANC_ALLELE' and 'DER_ALLELE' refer respectively to the ancestral site and the derived allele after the polarization. 
'CASE' and 'STATUS' inform us on whether the locus could be polarized. If it couldn't (due to missing data or an ambiguous allele configuration), 'STATUS' will display the string 'NOPOL'. Otherwise, the string will inform on which case scenario of allele distribution allowed us to polarize, along with the column 'CASE'. Note that the case scenario rarely matters to the user, in which case it's only a matter of whether the column 'STATUS' displays "NOPOL" (no polarization at this locus) or not (polarization at this locus).

--pop Name of the file containing the names of all the individuals in the target population included in the VCF, with one individual name per line.

--out1 Name of the file containing the names of all the individuals in the first outgroup population included in the VCF, with one individual name per line.

--out2 Name of the file containing the names of all the individuals in the second outgroup population included in the VCF, with one individual name per line.
