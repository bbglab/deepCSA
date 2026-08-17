process OMEGA_V2_RUN {
    tag "omega_v2"
    label 'cpu_medium'
    label 'process_high_memory'

    label 'bbgregressions'

    input:
    path(config)
    path(omega_assets)

    output:
    path("data.tsv")          , emit: data
    path("omega.tsv")         , emit: omega, optional: true
    path("omega.grouped.tsv") , emit: omega_grouped
    path "versions.yml"       , topic: versions

    script:
    def args = task.ext.args ?: ""
    """
    set -euo pipefail

    python ${omega_assets}/run.py \\
        --config ${config} \\
        --outfolder . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch data.tsv
    touch omega.tsv
    touch omega.grouped.tsv

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
