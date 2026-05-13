process OMEGA_MULTITEST {
    tag "all"
    label 'process_low'

    label 'deepcsa_core'

    input:
    path(omegas)
    tuple path(samples_json), path(groups_json), path(all_groups_json)

    output:
    path("${omegas.baseName}.tsv")  , emit: corrected
    path "versions.yml"            , topic: versions

    script:
    def output_name = "${omegas.baseName}.tsv"
    def temp_name = "${omegas.baseName}.corrected.tsv" // avoid overwriting staged input
    """
    omega_multiple_testing.py \\
        --omegas-file ${omegas} \\
        --samples-json ${samples_json} \\
        --groups-json ${groups_json} \\
        --output ${temp_name}

    mv ${temp_name} ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${omegas.baseName}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
