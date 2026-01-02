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
    tuple val(meta), path("*with_subgenic.tsv")  , emit: panel_increased
    tuple val(meta), path("subgenic_names.json") , emit: new_regions_json
    path "versions.yml"                          , topic: versions


    script:
    def autoexons = (params.autoexons && exons.size() > 0 ) ? "--autoexons ${exons}" : ""
    def autodomains = params.autodomains ? "--autodomains ${domains}" : ""
    def custom_regions = params.subgenic_bedfile ? "--custom ${custom}" : ""
    def subgenic_regions_complement = params.subgenic_regions_complement ? "--subgenic-regions-complement" : ""
    """
    add_subgenicregions.py --panel_file ${panel} \\
                        ${autoexons} \\
                        ${autodomains} \\
                        ${custom_regions} \\
                        ${subgenic_regions_complement}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch panel.with_subgenic.tsv
    touch subgenic_names.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
