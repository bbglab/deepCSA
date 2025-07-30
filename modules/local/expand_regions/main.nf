process EXPAND_REGIONS {

    tag "$meta.id"
    label 'process_high'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(panel)
    path (bedfile)

    output:
    tuple val(meta), path("*with_hotspots.tsv") , emit: panel_increased
    tuple val(meta), path("hotspot_names.json") , emit: new_regions_json
    path "versions.yml"                         , topic: versions


    script:
    def expansion = task.ext.expansion ?: 0
    def custom_bedfile = task.ext.using_bedfile ? "--bedfile custom_regions_bedfiles.bed" : ""
    def autoexons = params.omega_autoexons ? "--autoexons" : ""
    """
    cat ${bedfile} >> custom_regions_bedfiles.bed
    add_hotspots.py --panel_file ${panel} \\
                    ${custom_bedfile} \\
                    --expand ${expansion} \\
                    ${autoexons}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
