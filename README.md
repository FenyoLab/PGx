# PGx
PGx supports proteogenomic intregration of mass spectrometry (MS) proteomics data with next generation sequencing (NGS) by mapping identified peptides onto their putative genomic coordinates.

## Installation
The ideal installation setup for the python script-set which constitutes PGx is as a set of executable scripts in the path of the shell. For this you will need to move the scripts to a repository on the `$PATH` of your shell and add a `#!` line to the top of the script pointing to your preferred python installation. For example, if your python is installed in `/usr/bin/python` you should add `#!/usr/bin/python` to the top of every script, then `chmod` (see explanation for the chmod command [here](https://en.wikipedia.org/wiki/Chmod)) the scripts so they are executable (at which point it is customary to remove the `.py` suffix from the filename).

## Proteome Directories
PGx expects custom proteomes to be located in directories, where each directory contains a file called `proteome.fasta` and a file called `proteome.bed`. 

The file called `proteome.fasta` represent the sequence space that was searched by e.g. Mascot or Sequest during the proteomics phase of the Proteogenomics study. It must contain nothing but the entry name in the header line of each protein entry (the line that starts with: `">"`) and the protein sequence should span exactly one line (historically FASTA files have been limited to 60, 70 or 80 characters per line which is an unnecessary formatting requirement for modern bioinformatic pipelines. 

The file called `proteome.bed` represents a mapping of the sequence space back to its genomic coordinates. In proteogenomic studies this file is usually a side-effect of the sequencing effort and is often exported as a BED file (the format utilized by the PGx framework).

A minimal example proteome is provided in this repository (oncoproteome). The custom proteome contains an unexpected peptide discovered to be expressed in cancer cells as part of the [CPTAC](http://proteomics.cancer.gov/) research program (it is unexpected due to its presence in what is normally an intronic region of the gene).

## Example Usage
Once the scripts are installed, custom proteomes can be indexed and then queried for the location of peptides. The locations can then be mapped into a derived BED file which can be visualized as a custom track on all major genome browsers.

The sequence of commands is: 
```bash
pgx_index oncoproteome
pgx_query oncoproteome peptides.txt > hits.txt
pgx_bed oncoproteome hits.txt > hits.bed
```
If the scripts have not been installed as executable commands, the sequence can be run as follows:
```bash
python pgx_index.py oncoproteome
python pgx_query.py oncoproteome peptides.txt > hits.txt
python pgx_bed.py oncoproteome hits.txt > hits.bed
```
If the second argument is missing, the scripts assume the input will be piped into `STDIN`, which means it is possible to generate one liners such as:
```bash
echo SPPDSPTDALMQLAK | pgx_query oncoproteome | pgx_bed oncoproteome
```

## Miscellaneous Scripts
The misc directory contains additional scripts which support the automated download of RefSeq as well as the ability to query peptides with a possible nsSNP (essentially a single amino acid change in the protein sequence).
