process FILTER_BATCH {
    tag "$meta.id"

    label 'process_high_memory'
    label 'time_low'

    label 'deepcsa_core'

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.cohort.filtered.tsv.gz")  , emit: cohort_maf
    path "versions.yml"                                , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def repetitive_variant = task.ext.repetitive_variant ? "--repetitive-variant-threshold ${task.ext.repetitive_variant}" : ""
    def germline_threshold = task.ext.germline_threshold ? "--somatic-vaf-boundary ${task.ext.germline_threshold}" : ""
    def proportion_samples_nrich = task.ext.prop_samples_nrich ? "--n-rich-cohort-proportion ${task.ext.prop_samples_nrich}" : ""
    """
    filter_cohort.py \\
        --maf-df-file ${maf} \\
        --sample-name ${prefix} \\
        ${repetitive_variant} \\
        ${germline_threshold} \\
        ${proportion_samples_nrich}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch shared_cohort.filtered.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
