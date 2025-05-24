process ONCODRIVE3D_PREPROCESSING {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta),  path(maf)
    tuple val(meta2), path(all_vep_output)

    output:
    tuple val(meta), path("*.mutations.raw_vep.tsv") , emit: vep_output4o3d
    path "versions.yml"                              , topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    oncodrive3d_preprocessing.py ${maf} ${all_vep_output} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mutations.raw_vep.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
