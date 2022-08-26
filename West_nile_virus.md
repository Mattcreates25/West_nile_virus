# __Phylogenetic Analysis of the 2020 West Nile Virus (WNV) Outbreak in Andalusia (Spain)__

## introduction
The Andalusian viral samples belonged to lineage 1 and were relatively similar to those of previous outbreaks which occurred in the Mediterranean region

## workspace preparations 
The first thing we did was create a working directory where we will save our files. we initiated a git repo this way we could share our work on GitHub
 

```bash
mkdir WNV
cd WNV
git init
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
conda install snakemake
```
## Sequence retrieval 
download the sequences with wget
```bash
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953897.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953898.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953895.1?download=true
wget https://www.ebi.ac.uk/ena/browser/api/fasta/OU953896.1?download=true
```

combine the four sequences into a file that contains all the sequences then combined them with the ```WWR_sequences.fasta``` 
```bash
cat OU* >> four_sequences.fasta | cat WWR_sequences.fasta >> all_sequences.fasta
```

the reference sequence was downloaded from NCBI and saved into a file called ```refseq.fasta```
[link][1]

[1][https://www.ncbi.nlm.nih.gov/nuccore/NC_009942.1]

 the fasta sequences for the worldwide representative set of WNVs sequences were obtained using batch entrez and improted into a single file called 
 WWR_sequences.fasta
 accession numbers were retrieved from the ```viruses-13-00836-s001.zip``` which can be found in the ```Table S1.R3.xlsx.``` for retrieving purposes,
 they were saved into a text file called ```WWrep_accession.txt``` this text file was then uploaded into batch entrez

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

The entire analysis can be run using snakemake workflow management. ```Snakemake``` breaks a workflow into a set of rules that are specified in a file called Snakefile. Each rule takes a number of input files, specifies a few parameters, and produces output files

Parse delimited fields from FASTA sequence names into a TSV and FASTA file using a ```Snakefile```  input in this code creates a metadata file
alongside the fasta file and saved into a folder called results which will be used in subsequent steps
to do this we created a Snakemake file and used the augur parse command

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

create a ```python``` script to split the accessions from the description in the metadata and mark the script as an executable file

```bash
touch parse.py
chmod +x sample-script.py
```

convert the csv to tsv for augur filter
```bash
cat results/newmeta.csv | sed 's/,/\t/g' > results/newmeta.tsv
```

perform nucleotide base count  using ```augur index```
```bash
augur index -s data/all_sequences.fasta -o results/sequence_index.tsv
```

We then provide the sequence index as an input to augur filter commands to speed up filtering on sequence-specific attributes.
__encountered a problem here cause of a conflict in metadata__

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

for alignment ```mafft``` is required. MAFFT (Multiple Alignment using Fast Fourier Transform) is high-speeded multiple sequence alignment program.
```bash
sudo apt install mafft iqtree raxml fasttree vcftools
```

perform the MSA with Augur
```bash
align -s data/all_sequences.fasta -o results/all_alignment.fasta --method mafft --reference-sequence data/refseq.fasta --fill-gaps

```

generate tree
```bash
augur tree -a results/all_alignment.fasta -o results/alignment.nwk --method iqtree --substitution-model GTR -o alignment.nwk --tree-builder-args="-ninit 2 -n 2 -me 0.05"

```


# refine

```bash
augur refine -a results/all_alignment.fasta -t alignment.nwk --metadata results/all_metadata.tsv --timetree --output-tree results/refined_alignment.nwk --output-node-data results/branches.json

```

# export 
to export the tree due to metadata constraints as of now I used iTOL tree

![1 image](https://user-images.githubusercontent.com/97890823/186853004-83f14f42-99ef-47fd-abf5-9c56a327dcc4.png)
