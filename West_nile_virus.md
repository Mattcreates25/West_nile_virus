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

##combine the fasta files into one fasta file
cat OU* >> WNV_combi.fasta

### download the reference sequence from NCBI
[link][1]

[1][https://www.ncbi.nlm.nih.gov/nuccore/NC_009942.1]

## the fasta sequences for the world-wide representative set of WNVs sequences were obtained using batch entrez and improted into a single file

## for alignment mafft is required 
```bash
sudo apt install mafft iqtree raxml fasttree vcftools
```
## install figtree to visualize tree data
```
sudo apt-get -y install figtree
```
#perform the MSA with Augur
```bash
augur align [-h] --sequences FASTA [FASTA ...] [--output OUTPUT]
                   [--nthreads NTHREADS] [--method {mafft}]
                   [--reference-name NAME] [--reference-sequence PATH]
                   [--remove-reference] [--fill-gaps]
                   [--existing-alignment FASTA] [--debug]
```

```bash
augur align -s table_sequences.fasta WNV_combi.fasta  --method mafft --fill-gaps --reference-sequence NC_009942.1.fasta -o new_alignment.fasta

```

# generate tree
```bash
usage: augur tree [-h] --alignment ALIGNMENT
                  [--method {fasttree,raxml,iqtree}] [--output OUTPUT]
                  [--substitution-model SUBSTITUTION_MODEL]
                  [--nthreads NTHREADS] [--vcf-reference VCF_REFERENCE]
                  [--exclude-sites EXCLUDE_SITES]
                  [--tree-builder-args TREE_BUILDER_ARGS]
                  [--override-default-args]
```

```bash
augur tree -a new_alignment.fasta --method iqtree --substitution-model GTR -o alignment.nwk --tree-builder-args="-ninit 2 -n 2 -me 0.05"

```


