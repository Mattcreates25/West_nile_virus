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

##combined all the sequences ie the 4 reference sequences and the world-wide reperesentative set of WNVs
cat WNV_combi.fasta table_sequences.fasta >>all_sequences.fasta

## for alignment mafft is required 
```bash
sudo apt install mafft iqtree raxml fasttree vcftools
```
## install figtree to visualize tree data
```
sudo apt-get -y install figtree
```

## filter and subsample the sequences

usage: augur filter [-h] --metadata FILE [--sequences SEQUENCES]
                    [--sequence-index SEQUENCE_INDEX]
                    [--metadata-chunk-size METADATA_CHUNK_SIZE]
                    [--metadata-id-columns METADATA_ID_COLUMNS [METADATA_ID_COLUMNS ...]]
                    [--query QUERY] [--min-date MIN_DATE]
                    [--max-date MAX_DATE]
                    [--exclude-ambiguous-dates-by {any,day,month,year}]
                    [--exclude EXCLUDE [EXCLUDE ...]]
                    [--exclude-where EXCLUDE_WHERE [EXCLUDE_WHERE ...]]
                    [--exclude-all] [--include INCLUDE [INCLUDE ...]]
                    [--include-where INCLUDE_WHERE [INCLUDE_WHERE ...]]
                    [--min-length MIN_LENGTH] [--non-nucleotide]
                    [--group-by GROUP_BY [GROUP_BY ...]]
                    [--sequences-per-group SEQUENCES_PER_GROUP | --subsample-max-sequences SUBSAMPLE_MAX_SEQUENCES]
                    [--probabilistic-sampling | --no-probabilistic-sampling]
                    [--priority PRIORITY] [--subsample-seed SUBSAMPLE_SEED]
                    [--output OUTPUT] [--output-metadata OUTPUT_METADATA]
                    [--output-strains OUTPUT_STRAINS]
                    [--output-log OUTPUT_LOG]

#perform the MSA with Augur
```bash
augur align [-h] --sequences FASTA [FASTA ...] [--output OUTPUT]
                   [--nthreads NTHREADS] [--method {mafft}]
                   [--reference-name NAME] [--reference-sequence PATH]
                   [--remove-reference] [--fill-gaps]
                   [--existing-alignment FASTA] [--debug]
```

```bash
augur align -s all_sequences.fasta --method mafft --fill-gaps --reference-sequence NC_009942.1.fasta -o new_alignment.fasta


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


# refine
```bash
usage: augur refine [-h] [--alignment ALIGNMENT] --tree TREE [--metadata FILE]
                    [--output-tree OUTPUT_TREE]
                    [--output-node-data OUTPUT_NODE_DATA] [--use-fft]
                    [--timetree] [--coalescent COALESCENT]
                    [--gen-per-year GEN_PER_YEAR] [--clock-rate CLOCK_RATE]
                    [--clock-std-dev CLOCK_STD_DEV] [--root ROOT [ROOT ...]]
                    [--keep-root] [--covariance] [--no-covariance]
                    [--keep-polytomies] [--precision {0,1,2,3}]
                    [--date-format DATE_FORMAT] [--date-confidence]
                    [--date-inference {joint,marginal}]
                    [--branch-length-inference {auto,joint,marginal,input}]
                    [--clock-filter-iqd CLOCK_FILTER_IQD]
                    [--vcf-reference VCF_REFERENCE]
                    [--year-bounds YEAR_BOUNDS [YEAR_BOUNDS ...]]
                    [--divergence-units {mutations,mutations-per-site}]
                    [--seed SEED]
```

```bash

```
