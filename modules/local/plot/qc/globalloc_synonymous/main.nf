process PLOT_OMEGASYN_QC {

    tag "all"
    label 'process_low'

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"

    input:
    path(all_obs_syn_muts)
    path(all_syn_muts_gloc)
    tuple path (samples_json), path (groups_json), path (all_groups_json)

    output:
    path("**.pdf")  , emit: plots_pdf
    path("**.png")  , emit: plots_png
    path("**.tsv")  , emit: tables

    path "versions.yml"              , topic: versions


    script:
    def prefix = task.ext.prefix ?: "omega_qc"

    """
    mkdir ${prefix}.plots
    ls -1 ${all_obs_syn_muts} > observed.txt
    ls -1 ${all_syn_muts_gloc} > estimated.txt
    omega_syn_qc.py \\
                    --observed-syn observed.txt \\
                    --estimated-syn estimated.txt \\
                    --output-prefix ${prefix}.plots/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
