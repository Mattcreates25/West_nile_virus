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
        sequences = rules.parse.output.sequences,
	reference = "data/refseq.fasta"
    params:
        method = "mafft"
    output:
        alignment="results/all_alignment.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
	    --reference-sequence {input.reference} \
            --method {params.method} \
            --output {output.alignment} \
	    --fill-gaps
        """
rule tree:
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method iqtree \
            --substitution-model GTR \
            --tree-builder-args="-ninit 2 -n 2 -me 0.05"
        """

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = "results/newmeta.csv"
    output:
        tree = "results/new_tree.nwk",
        node_data = "results/branches.json"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """

rule export:
    input:
        tree  = rules.refine.output.tree,
	node_data = rules.refine.output.node_data
    output:
        tree = "results/branches_na.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
	    --node-data {input.node_data} \
            --output {output.tree} \
	    --skip-validation
        """

