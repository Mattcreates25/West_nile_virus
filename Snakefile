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

rule index:
    input:
        sequences = rules.parse.output.sequences
    output:
        indexed = "results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output-sequences {output.indexed}
        """


rule align:
    input:
        sequences = rules.parse.output.sequences
	reference = data/req_seq.fasta
    params:
        method = mafft
    output:
        alignment="results/all_alignment.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
	    --reference-sequence {input.reference}
            --method {params.method} \
            --output {output.alignment} \
	    --fill-gaps
        """

rule tree:
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree.nwk"
    params:
	methods = iqtree
	substmodel = GTR
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.methods} \
	    --substitution-model {params.substmodel} \
	    --tree-builder-args="-ninit 2 -n 2 -me 0.05"
        """
