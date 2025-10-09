process EXPAND_REGIONS {

    tag "$meta.id"
    label 'process_high'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    tuple val(meta), path(panel)
    path (domains)
    path (exons)
    path (custom)

    output:
    tuple val(meta), path("*with_hotspots.tsv") , emit: panel_increased
    tuple val(meta), path("hotspot_names.json") , emit: new_regions_json
    path "versions.yml"                         , topic: versions


    script:
    // def expansion = task.ext.expansion ?: 0
    def autoexons = params.omega_autoexons ? "--autoexons ${exons}" : ""
    def autodomains = params.omega_autodomains ? "--autodomains ${domains}" : ""
    def custom_regions = params.omega_subgenic_bedfile ? "--custom ${custom}" : ""
    """
    add_hotspots.py --panel_file ${panel} \\
                        ${autoexons} \\
                        ${autodomains} \\
                        ${custom_regions}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch panel.with_hotspots.tsv
    touch hotspot_names.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
