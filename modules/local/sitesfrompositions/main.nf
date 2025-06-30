process SITESFROMPOSITIONS {

    tag "${meta.id}"

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"


    input:
    tuple val(meta), path(depths)

    output:
    tuple val(meta), path("*.sites4VEP.tsv")  , emit: annotated_panel_reg
    path  "versions.yml"                      , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assembly = task.ext.assembly ?: "hg38"
    // TODO
    // see if there is a better way to filter out chromosomes
    // that are not the canonical ones right now doing it by
    // filtering the size of the chromosome name to be smaller of equal to 2
    """
    cat <(printf "CHROM\\tPOS\\n") <( zcat ${depths} | cut -f1,2  | sed 's/^chr//g' | awk 'length(\$1) <= 2' ) > captured_positions.tsv;
    sites_table_from_positions.py \\
                    captured_positions.tsv \\
                    ${assembly} \\
                    captured_positions.sites4VEP.tmp.tsv;

    awk '{print "chr"\$0}' captured_positions.sites4VEP.tmp.tsv > captured_positions.sites4VEP.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch captured_positions.sites4VEP.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
