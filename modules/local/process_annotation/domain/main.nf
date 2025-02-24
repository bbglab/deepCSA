process DOMAIN_ANNOTATION {

    tag "${meta.id}"

    label 'cpu_low'
    label 'time_low'
    label 'process_high_memory'

    container "docker.io/ferriolcalvet/bgreference"

    input:
    tuple val(meta) , path(panel_annotated)
    path (domains_file)

    output:
    path("*.domains.bed4.bed")  , emit: domains_bed
    path  "versions.yml"        , topic: versions

    when:
    task.ext.when == null || task.ext.when

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
