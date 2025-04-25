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

    when:
    task.ext.when == null || task.ext.when

    script:
    def expansion = task.ext.expansion ?: 0
    def bedfile_to_use = task.ext.using_bedfile ? "custom_regions_bedfiles.bed" : "None"
    def autoexons = params.omega_autoexons ? "1" : "0"
    """
    cat ${bedfile} >> custom_regions_bedfiles.bed
    add_hotspots.py ${panel} \\
                    ${bedfile_to_use} \\
                    ${expansion} \\
                    ${autoexons}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
