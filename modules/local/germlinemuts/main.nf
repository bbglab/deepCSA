process GERMLINE_MUTATIONS {

    tag "$meta.id"
    label 'deepcsa_core'

    input:
    tuple val(meta),  path(clean_maf)
    tuple val(meta2), path(somatic_maf)

    output:
    tuple val(meta), path("*.germline.mutations.tsv"), emit: germline_mutations
    tuple val(meta), path("*pathogenic_snps*.tsv")   , emit: pathogenic_snps
    tuple val(meta), path("*.ancestry_inference.tsv"), emit: ancestry
    tuple val(meta), path("*.pdf")                   , emit: plots
    path "versions.yml"                              , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def gnomad_af_threshold = task.ext.gnomad_af_threshold ? "--gnomad-af-threshold ${task.ext.gnomad_af_threshold}" : ""
    """
    summarize_germline_snps.py \\
        --clean_maf ${clean_maf} \\
        --somatic_maf ${somatic_maf} \\
        --output_prefix ${prefix} \\
        ${gnomad_af_threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.germline.mutations.tsv
    touch ${prefix}.pathogenic_snps.tsv
    touch ${prefix}.pathogenic_snps_summary.tsv
    touch ${prefix}.ancestry_inference.tsv
    touch ${prefix}.germline_snps_summary.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
