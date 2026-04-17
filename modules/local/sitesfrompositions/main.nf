process SITESFROMPOSITIONS {

    tag "${meta.id}"

    label 'deepcsa_core'

    input:
    tuple val(meta), path(depths)

    output:
    tuple val(meta), path("*.sites4VEP.chunk*.tsv")  , emit: annotated_panel_reg
    path  "versions.yml"                             , topic: versions


    script:
    def assembly = task.ext.assembly ?: "hg38"
    def chunk_size = task.ext.chunk_size ?: 0

    // TODO
    // see if there is a better way to filter out chromosomes
    // that are not the canonical ones right now doing it by
    // filtering the size of the chromosome name to be smaller of equal to 2

    """
    cat <(printf "CHROM\\tPOS\\n") <( zcat ${depths} | cut -f1,2  | sed 's/^chr//g' | awk 'length(\$1) <= 2' ) > captured_positions.tsv;
    sites_table_from_positions.py \\
        --input-positions captured_positions.tsv \\
        --genome-assembly ${assembly} \\
        --output-prefix captured_positions.sites4VEP \\
        --chunk-size ${chunk_size};

    rm captured_positions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch captured_positions.sites4VEP.chunk0.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
