process ADAPT_PANEL_REFCDS {

    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'

    label 'deepcsa_core'

    input:
    tuple val(meta) , path(biomart)
    tuple val(meta2), path(bedfile)

    output:
    tuple val(meta), path("custom_filtered_biomart.tsv"), emit: filtered_biomart
    tuple val(meta), path("splice_sites.tsv")           , emit: splice_sites
    path "versions.yml"                                 , topic: versions

    script:
    """
    dNdScv_panel_prep.py --bed ${bedfile} \\
                        --genes ${biomart} \\
                        --output custom_filtered_biomart.tsv \\
                        --verbose \\
                        -s splice_sites.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dNdScv_panel_prep : \$(dNdScv_panel_prep.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch custom_filtered_biomart.tsv  
    touch splice_sites.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dNdScv_panel_prep : \$(dNdScv_panel_prep.py --version)
    END_VERSIONS
    """
}