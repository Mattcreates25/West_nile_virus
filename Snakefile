rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = "data/metadata.csv"
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --timetree \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """
