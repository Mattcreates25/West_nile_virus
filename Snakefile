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

rule align:
    input:
        seq="data/all_sequences.fasta"
    params:
        nthreads = 2
    output:
        aln="results/alignment.fasta"
    shell:
        '''
        augur align --sequences {input.seq} --nthreads {params.nthreads} --ouput {output.aln}
        '''
rule tree:
    input:
        aln="results/all_alignment.fasta"
    output:
        tree="results/tree.nwk"
    shell:
        '''
        augur tree --alignment {input.aln} --ouput {output.tree}
        '''
