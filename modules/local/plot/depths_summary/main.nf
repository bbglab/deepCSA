process PLOT_DEPTHS {

    tag "$meta.id"

    label 'deepcsa_core'

    input:
    tuple val(meta) , path(depth)
    tuple val(meta2), path(panel)
    path(region_file)

    output:
    tuple val(meta), path("*.pdf")             , optional : true, emit: plots
    tuple val(meta), path("*.avgdepth_per_sample.tsv")          , emit: average_per_sample
    tuple val(meta), path("*.avgdepth_per_gene.tsv")            , emit: average_per_gene
    tuple val(meta), path("*.depth_per_gene_per_sample.tsv")    , emit: average_per_gene_sample
    tuple val(meta), path("*.avgdepth_per_region.tsv")          , optional : true, emit: average_per_region
    tuple val(meta), path("*.depth_per_region_per_sample.tsv")  , optional : true, emit: average_per_region_sample
    tuple val(meta), path("*.genes_per_region.tsv")             , optional : true, emit: genes_per_region
    tuple val(meta), path("*depth*.tsv")       , optional : true, emit: depths
    path  "versions.yml"                                        , topic: versions



    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    def plot_within_gene = task.ext.withingene ? "True" : "False"
    // Only pass the region file when a real one is provided (a placeholder is used when none is set)
    def region_arg = params.chromosome_region_file ? "--region_file ${region_file}" : ""
    """
    plot_depths.py \\
                --sample_name ${prefix} \\
                --depth_file ${depth} \\
                --panel_bed6_file ${panel} \\
                --panel_name ${panel_version} \\
                --plot_within_gene ${plot_within_gene} \\
                ${region_arg};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    touch ${prefix}.${panel_version}.depths_info.tsv
    touch ${prefix}.avgdepth_per_region.tsv
    touch ${prefix}.depth_per_region_per_sample.tsv
    touch ${prefix}.genes_per_region.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
