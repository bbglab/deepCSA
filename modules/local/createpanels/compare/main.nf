process COMPARE_TRINUCLEOTIDE_PROPORTIONS_PANELS {
    tag "panels"

    label 'deepcsa_core'


    input:
    path (consensus_panels)
    path (wgs_trinucleotides)


    output:
    path ("*.png")          , emit: plots
    path  "versions.yml"    , topic: versions


    script:
    """
    qc_trinucleotide_proportions.py --wgs-trinucleotide ${wgs_trinucleotides}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.mapping.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
