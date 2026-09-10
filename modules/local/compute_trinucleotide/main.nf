process COMPUTE_TRINUCLEOTIDE {

    tag "$meta.id"
    label 'cpu_low'

    label 'deepcsa_core'

    input:
    tuple val(meta), path(depths)
    path (wgs_trinuc)

    output:
    tuple val(meta), path("*.trinucleotides.tsv.gz")                    , emit: trinucleotides
    tuple val(meta), path("*.png")                  , optional : true   ,  emit: proportions_plot
    path "versions.yml"                                                 , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    mutprof_2compute_trinucleotide.py \\
                    --depths_file ${depths} \\
                    --sample_name ${prefix} \\
                    --ref_wgs_trinucleotides ${wgs_trinuc} \\
                    ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """


    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.profile.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
