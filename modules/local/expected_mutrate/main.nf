process EXP_MUTRATE {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/axelrosendahlhuber/expected_mutrate:latest'

    input:
    tuple val(meta), path(regions)
    tuple val(meta2), path(mutations)
    tuple val(meta3), path(depths)
    tuple val(meta4), path(annotated_panel)
    tuple val(meta5), path(annotated_bed_file)
    path (features_table)


    output:
    tuple val(meta), path("**.png")     , emit: plots
    tuple val(meta), path("**.tsv")     , emit: stats
    tuple val(meta), path("**.rds")     , emit: rds_file
    tuple val(meta), path("**.rda")     , emit: rda_file

    path "versions.yml"                 , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metadata_file = task.ext.metadata_file ? "${features_table}": ""
    """
    mkdir expected_mutrate
    mutrisk_deepCSA.R \\
            ${regions} \\
            ${mutations} \\
            ${depths} \\
            ${annotated_panel} \\
            expected_mutrate \\
            ${annotated_bed_file} \\
            ${metadata_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
    END_VERSIONS
    """
}
