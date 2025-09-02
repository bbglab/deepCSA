process PLOT_METRICS_QC {

    tag "$meta.id"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    path(all_mutdensities)
    tuple val(meta2), path(panel_file)

    path (samples_json)

    output:
    tuple val(meta), path("**.pdf")  , emit: plots
    tuple val(meta), path("**.csv")  , emit: tables

    path "versions.yml"              , topic: versions


    script:
    def prefix = task.ext.prefix ?: "all_samples"

    """
    mkdir ${prefix}.plots
    mutation_densities_qc.py \\
                    --input-file ${all_mutdensities} \\
                    --output-dir ${prefix}.plots \\
                    --panel ${panel_file} \\
                    --samples ${samples_json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
