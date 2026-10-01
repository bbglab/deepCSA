process DNA_2_PROTEIN_MAPPING {
    tag "$meta.id"

    label 'deepcsa_core'


    input:
    tuple val(meta) , path(mutations_file)
    tuple val(meta2), path(panel_file)
    tuple val(meta3), path(all_samples_depths)
    path(gff3_file)


    output:
    tuple val(meta2), path("depths_per_position_exon_gene.tsv") , emit: depths_exons_positions
    tuple val(meta2), path("panel_exons.bed4.bed")              , emit: panel_exons_bed
    tuple val(meta2), path("*.pdf")                             , emit: covered_genes_pdf
    tuple val(meta2), path("coverage*.tsv")                     , emit: covered_genes_tsv
    path  "versions.yml"                                        , topic: versions


    script:
    def ensembl_release = "--ensembl-release \"${task.ext.ensembl_release}\""
    def ensembl_species = "--ensembl-species \"${task.ext.ensembl_species}\""
    def ensembl_genome  = "--ensembl-genome \"${task.ext.ensembl_genome}\""
    def gff3_local = task.ext.gff3_provided ? "--gff3-file \"${gff3_file}\"" : ""
    """
    cut -f 1,2,6 ${panel_file} | uniq > ${meta2.id}.panel.unique.tsv
    panels_computedna2protein.py \\
                --mutations-file ${mutations_file} \\
                --consensus-file ${meta2.id}.panel.unique.tsv \\
                --depths-file ${all_samples_depths} \\
                ${ensembl_release} \\
                ${ensembl_species} \\
                ${ensembl_genome} \\
                ${gff3_local}


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
