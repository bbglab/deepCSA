process COMPARE_TRINUCLEOTIDE_PROPORTIONS_PANELS {
    tag "1"
    label 'process_single'

    label 'deepcsa_core'


    input:
    path (consensus_panels)
    path (wgs_trinucleotides)


    output:
    path ("*.png")          , emit: plots
    path  "versions.yml"    , topic: versions


    script:
    """
    compare_trinucleotide_proportions.py --wgs-trinucleotide ${wgs_trinucleotides}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
