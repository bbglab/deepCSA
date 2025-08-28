process DOMAIN_ANNOTATION {

    tag "${meta.id}"

    label 'cpu_low'
    label 'time_low'
    label 'process_high_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta) , path(panel_annotated)
    path (domains_file)

    output:
    path("*.domains.bed4.bed")  , emit: domains_bed
    path("domains_info.tsv")    , emit: domains_tsv
    path  "versions.yml"        , topic: versions


    script:
    """
    panels_dna2domain.py \\
                    --consensus-panel-rich ${panel_annotated} \\
                    --domains ${domains_file} \\
                    --output seq_panel.domains.bed4.bed;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch consensus.domains.bed4.bed;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
