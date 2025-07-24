process COMPUTE_CONTAMINATION {

    tag "$meta.id"
    label 'process_high'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(maf)
    tuple val(meta2), path(somatic_maf)


    output:
    tuple val(meta), path("*.tsv")  ,                   emit: contamination_results
    tuple val(meta2), path("*.pdf") ,   optional:true,  emit: contamination_plots
    path "versions.yml" ,                               topic: versions

    script:
    """
    check_contamination.py \\
                        --maf_path ${maf} \\
                        --somatic_maf ${somatic_maf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch contaminated_samples.tsv
    touch contaminated_samples.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
