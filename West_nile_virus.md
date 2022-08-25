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

## install auspice
```bash
conda install -c conda-forge nodejs
npm install --global auspice
```

## install snakemake
```bash
conda install snakemake
```

## download the sequences with wget
```bash
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953897.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953898.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953895.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953896.1?download=true
```

##combine all the sequences 
```bash
cat WWR_sequences.fasta WWR_sequences.fasta >>all_sequences.fasta
```

###  the reference sequence was downloaded from NCBI and saved into a file called refseq.fasta
[link][1]

[1][https://www.ncbi.nlm.nih.gov/nuccore/NC_009942.1]

## the fasta sequences for the world-wide representative set of WNVs sequences were obtained using batch entrez and improted into a single file called 
## WWR_sequences.fasta
## accession numbers were retrieved from the viruses-13-00836-s001.zip which can be found in the Table S1.R3.xlsx.for retrieving purposed they were saved into a 
## text file called WWrep_accession.txt
```bash
wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8148183/bin/viruses-13-00836-s001.zip
unzip viruses-13-00836-s--1.zip
```
## create a txt file and paste all the accession numbers there. this file will be input into batch entrez
```bash
touch WWrep_accession.txt
```


## create a directory for all sequence data
```bash
mkdir data
```

#phylogenetic analysis

##Parse delimited fields from FASTA sequence names into a TSV and FASTA file using a Snakefile  input this code this creates a metadata file alongside the fasta file
## and saved into a folder called results which will be used in subsequent steps
```bash
rule parse:
    input:
        sequences = "data/all_sequences.fasta"
    output:
        sequences = "results/all_sequencesP.fasta",
        metadata = "results/all_metadata.tsv"
    params:
        fields = "site accession strain desciption"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --fields {params.fields} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
```



## perform nucleotide base count  using augur index
```bash
augur index -s data/all_sequences.fasta -o results/sequence_index.tsv
```

## We then provide the sequence index as an input to augur filter commands to speed up filtering on sequence-specific attributes.

```bash
augur filter 
```

## for alignment mafft is required 
```bash
sudo apt install mafft iqtree raxml fasttree vcftools
```

#perform the MSA with Augur
```bash
align -s data/all_sequences.fasta -o results/all_alignment.fasta --method mafft --reference-sequence data/refseq.fasta --fill-gaps

```

# generate tree
```bash
augur tree -a results/all_alignment.fasta -o results/alignment.nwk --method iqtree --substitution-model GTR -o alignment.nwk --tree-builder-args="-ninit 2 -n 2 -me 0.05"

```


# refine

```bash
augur refine -a results/all_alignment.fasta -t alignment.nwk --metadata results/all_metadata.tsv --timetree --output-tree results/refined_alignment.nwk --output-node-data results/branches.json

```

#export 

