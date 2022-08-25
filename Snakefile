rule parse:
    input:
        sequences = "data/all_sequences.fasta"
    output:
        sequences = "results/all_sequencesP.fasta",
        metadata = "results/all_metadata.tsv"
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
