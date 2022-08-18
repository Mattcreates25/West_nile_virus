# Phylogenetic Analysis of the 2020 West Nile Virus (WNV) Outbreak in Andalusia (Spain)

## make a directory for your work and move into said directory while you are at it initialize a git repo

```bash
mk wdir WNV
cd WNV
git init
```
## Create a conda environment to install augur into and activate that environment go ahead and install augur and its dependencies into your environment. 
```bash
conda create -n forAugur
conda activate forAugur
conda install -c conda-forge -c bioconda augur
```

## download the sequences with wget
```
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953897.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953898.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953895.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953896.1?download=true
```
### download the reference sequence from NCBI
[link][1]

[1][https://www.ncbi.nlm.nih.gov/nuccore/NC_009942.1]

## for alignment mafft is required 
```bash
sudo apt install mafft iqtree raxml fasttree vcftools
```

#perform the MSA with Augur
```augur align [-h] --sequences FASTA [FASTA ...] [--output OUTPUT]
                   [--nthreads NTHREADS] [--method {mafft}]
                   [--reference-name NAME] [--reference-sequence PATH]
                   [--remove-reference] [--fill-gaps]
                   [--existing-alignment FASTA] [--debug]
```


