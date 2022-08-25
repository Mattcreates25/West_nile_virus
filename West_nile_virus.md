# __Phylogenetic Analysis of the 2020 West Nile Virus (WNV) Outbreak in Andalusia (Spain)__

##workspace preparations 
The first thing we did was create a working directory where we will save our files. we initiated a git  repo this way we could share our work on GitHub
 

```bash
mk wdir WNV
cd WNV
git init
```
### Create a conda environment
Augur is a bioinformatic tool for pylogenetic analysis. the collection of comands from the tool are designed to be used with a larger processing pipeline like
```snakemake``` Augur is composed of a series of modules and different workflows will use different parts of the pipeline.
 A selection of augur modules and different possible entry points are illustrated below.


![Augur](https://docs.nextstrain.org/projects/augur/en/stable/_images/augur_analysis_sketch.png)

For Augur to run we created a unique environment where we could install the package using miniconda 
###install augur
```bash
conda create -n forAugur
conda activate forAugur
conda install -c conda-forge -c bioconda augur
```
Nextstrains auspice is an open source interactive tool for vizualizing phylogenetic data
### install auspice
```bash
conda install -c conda-forge nodejs
npm install --global auspice
```
snakemake is the pipeline tool preferred by nextstrain in simplifying augur commands
### install snakemake
```bash
conda install snakemake
```
##Sequence retrieval 
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

 the fasta sequences for the world-wide representative set of WNVs sequences were obtained using batch entrez and improted into a single file called 
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

##phylogenetic analysis

Parse delimited fields from FASTA sequence names into a TSV and FASTA file using a ```Snakefile```  input in this code creates a metadata file
alongside the fasta file and saved into a folder called results which will be used in subsequent steps
to do this we created a Snakemake file  and used the augur parse command
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


perform nucleotide base count  using ```augur index```
```bash
augur index -s data/all_sequences.fasta -o results/sequence_index.tsv
```

We then provide the sequence index as an input to augur filter commands to speed up filtering on sequence-specific attributes.
__encountered a problem here cause of conflict in metadata__

```bash
augur filter --sequences data/all_sequences.fasta --metadata results/metadata.tsv --sequence-index results/sequence_index.tsv  --output filtered.fasta 
```

for alignment ```mafft``` is required. MAFFT (Multiple Alignment using Fast Fourier Transform) is a high speed multiple sequence alignment program.
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

#export 

