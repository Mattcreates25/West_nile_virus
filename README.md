# __Phylogenetic Analysis of the 2020 West Nile Virus (WNV) Outbreak in Andalusia (Spain)__

## introduction

West Nile virus (WNV), a member of the Flavivirus genus, is transmitted in an enzootic cycle involving birds as amplifying hosts and mosquitoes as vectors, which can ultimately be transmitted to mammals, considered dead-end hosts, causing disease outbreaks in horses and/or humans. Currently, the virus is considered a recurrent zoonosis with a wide geographic distribution. Phylogenetically, WNV is classified into eight lineages.The Andalusian viral samples belonged to lineage 1 and were relatively similar to those of previous outbreaks which occurred in the Mediterranean region.
A phylogenetic analysis was performed on the obtained consensus genomes in the context of a world-wide representative set of WNVs. 

## Our objectives
1. Download sequences
2. Phylogenetics Analysis
3. Document your work in GitHub
4. Prepare a presentation for the same

## workspace preparations 
The first thing we did was create a working directory where we will save our files. we initiated a git repo this way we could share our work on GitHub
 

```bash
mkdir WNV
cd WNV
git init
```
you can copy this github repo using ```git clone``` in your terminal

```bash
git clone git@github.com:Mattcreates25/West_nile_virus.git
```

### Create a conda environment
Augur is a bioinformatics tool for phylogenetic analysis. the collection of commands from the tool are designed to be used with a larger processing pipeline like
```snakemake``` Augur is composed of a series of modules and different workflows will use different parts of the pipeline.
 A selection of augur modules and different possible entry points are illustrated below.


![Augur](https://docs.nextstrain.org/projects/augur/en/stable/_images/augur_analysis_sketch.png)

For Augur to run we created a unique environment where we could install the package using miniconda 

### install augur
```bash
conda create -n forAugur
conda activate forAugur
conda install -c conda-forge -c bioconda augur
```
Nextstrain’s auspice is an open-source interactive tool for visualizing phylogenetic data

### install auspice
```bash
conda install -c conda-forge nodejs
npm install --global auspice
```
snakemake is the pipeline tool preferred by Nextstrain in simplifying augur commands

### install snakemake
```bash
conda install snakemake -c <channel>
```
## Sequence retrieval 
download the sequences with wget
```bash
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953897.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953898.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953895.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953896.1?download=true
```

The reference sequence was downloaded from NCBI and saved into a file called ```refseq.fasta```

[1][https://www.ncbi.nlm.nih.gov/nuccore/NC_009942.1]

 The fasta sequences for the worldwide representative set of WNVs sequences were obtained using batch entrez and improted into a single file called 
 WWR_sequences.fasta
 accession numbers were retrieved from the ```viruses-13-00836-s001.zip``` which can be found in the ```Table S1.R3.xlsx.``` for retrieving purposes,
 they were saved into a text file called ```WWrep_accession.txt``` this text file was then uploaded into batch entrez.
 
combine the four sequences into a file that contains all the sequences then combined them with the ```WWR_sequences.fasta``` 
```bash
cat OU* >> four_sequences.fasta | cat WWR_sequences.fasta >> all_sequences.fasta
```

```bash
wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8148183/bin/viruses-13-00836-s001.zip
unzip viruses-13-00836-s--1.zip
```
## create a txt file and paste all the accession numbers there. this file will be input into batch entrez
```bash
touch WWrep_accession.txt
```


create a directory for all sequence data and move all the sequence data into that directory
```bash
mkdir data 
mv *.fasta data
```

## phylogenetic analysis

The entire analysis can be run using snakemake workflow management. ```Snakemake``` breaks a workflow into a set of rules that are specified
in a file called Snakefile. Each rule takes a number of input files, specifies a few parameters, and produces output files

Parse delimited fields from FASTA sequence names into a TSV and FASTA file using a ```Snakefile```  input in this code creates a metadata file
alongside the fasta file and saved into a folder called results which will be used in subsequent steps
to do this we created a Snakemake file and used the augur parse command

```bash
rule parse:
    input:
        sequences = "data/all_sequences.fasta"
    output:
        sequences = "results/all_sequencesP.fasta",
        metadata = "results/all_metadata.csv"
    params:
        fields = "strain"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --fields {params.fields} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata}
        """

```

we Created a ```python``` script to split the accessions from the description in the metadata and mark the script as an executable file

```py
# import the panda library
import pandas as pd

#read in the metadata csv
df = pd.read_csv('results/all_metadata.csv')

#check the format
df.head(15)

#create a variable with split data
new_metadata = df['strain'].str.split('.', n=1, expand=True)

#rename the columns
new_meta = new_metadata.rename(columns={0:'strain',1:'name'})

#export as a csv
new_meta.to_csv('newmeta.csv',index=False)

```

__marking the parse script as an executable file__
```bash
touch parse.py
chmod +x sample-script.py
```

perform nucleotide base count  using ```augur index```
```bash
augur index -s data/all_sequences.fasta -o results/sequence_index.tsv
```

We then provide the sequence index as an input to augur filter commands to help in the filtering of sequence-specific attributes.
__encountered a problem here cause of a conflict in metadata opted to skip the filter option and worked with the entire data set__

## error message
```bash
augur filter --sequences data/all_sequences.fasta --metadata results/all_metadata.csv --sequence-index results/sequence_index.tsv  --output filtered.fasta 
```
```bash
Traceback (most recent call last):
  File "/home/icipe/miniconda3/envs/forAugur/lib/python3.10/site-packages/augur/__init__.py", line 65, in run
    return args.__command__.run(args)
  File "/home/icipe/miniconda3/envs/forAugur/lib/python3.10/site-packages/augur/filter.py", line 1400, in run
    metadata_reader = read_metadata(
  File "/home/icipe/miniconda3/envs/forAugur/lib/python3.10/site-packages/augur/io/metadata.py", line 76, in read_metadata
    raise Exception(f"None of the possible id columns ({id_columns!r}) were found in the metadata's columns {tuple(chunk.columns)!r}")
Exception: None of the possible id columns (['strain', 'name']) were found in the metadata's columns ('s', 'rain')


An error occurred (see above) that has not been properly handled by Augur.
To report this, please open a new issue including the original command and the error above:
    <https://github.com/nextstrain/augur/issues/new/choose>

```

for alignment ```mafft``` is required. MAFFT (Multiple Alignment using Fast Fourier Transform) and it is a high-speed multiple sequence alignment program.
```bash
sudo apt install mafft iqtree raxml fasttree vcftools
```

## perform the MSA with Augur
```bash
augur align -s data/all_sequences.fasta -o results/all_alignment.fasta --method mafft --reference-sequence data/refseq.fasta --fill-gaps

```

## generate tree
The phylogenetic tree was recovered by maximum likelihood, using a general time reversible model ```Augur version 14 +``` introduced --tree-builder-args which allows the user to create a tree with bootstrap values.

```bash
augur tree -a results/all_alignment.fasta -o results/tree.nwk --method iqtree --substitution-model GTR -o results/tree.nwk --tree-builder-args="-ninit 2 -n 2 -me 0.05"
```


## refine
Branching date estimation was carried out with the least square dating (LSD2) method which is the default.

```bash
augur refine --tree results/tree.nwk -a results/all_alignment.fasta --metadata results/newmeta.csv --output-tree results/new_tree.nwk --output-node-data results/branches.json --keep-root
```

## export
```bash
augur export v2 -t results/new_tree.nwk --node-data results/branches.json --output results/tree_auspice.json
``` 

# view the tree 
To view the tree we will have to write a ```narrative.md``` file that works alongside the dataset directory. This can be found in  ```narratives``` 

```
---
title: west_Nile_narrative
authors: "Mark Njama"
authorLinks: "https://github.com/Mattcreates25"
affiliations: "icipe"
date: "August 2022"
dataset: "http://localhost:4000/results/branches/na?d=tree"
abstract: "This narrative will take us to auspice for vizualization."
---

```

after this you can view your tree by running this line in bash this creates a localhost with your tree(click on the ```localhost 4000``` link)

```bash
auspice view --datasetDir results/ --narrativeDir narrative/
 
```

![paperstree](https://user-images.githubusercontent.com/97890823/187367745-caa5f985-3f3e-4bb5-8da8-675517646ea7.jpg)

![personaltree](https://user-images.githubusercontent.com/97890823/187367580-2477f803-d3d9-4206-961e-4b3c2937a40d.png)




Sequences of the Spanish 2020 WNV outbreak ```OU``` the closest relatives from previous outbreaks in Italy and the sequence JF719069 from a lethal equine case in Andalusia (Spain) in 2010. Other Spanish outbreaks were: JF707789, from a mosquito in Huelva, FJ766331 and FJ766332 from a golden eagle in Toledo. Other related outbreaks from the Mediterranean region (Cyprus MF797870), or adjacent locations (United Arab Emirates KU588135 and Russia MN149538) are also included.

## future
use UFboot to  generate bootstrap values
calculate non-synonymous to synonymous ratios along the viral genomes using the ```KaKs_Calculator```
