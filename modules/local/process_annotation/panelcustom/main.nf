process CUSTOM_ANNOTATION_PROCESSING {

    tag "${meta.id}"

    label 'cpu_low'
    label 'time_low'
    label 'process_high_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(panel_annotated)
    path custom_regions

    output:
    tuple val(meta), path("*.custom.tsv"), emit: custom_panel_annotation
    tuple val(meta), path("added_regions.tsv"), emit: added_regions
    path "versions.yml", topic: versions

    script:
    def simple = task.ext.simple ? "True" : "False"
    // TODO
    // Document this custom_regions has to be a TSV file with the following columns:
    // chromosome  start   end gene_name    impactful_mutations [neutral_impact] [new_impact]
    // chromosome start and end indicate the region that is being customized
    // gene_name           : is the name of the region that is being added, make sure that it does not coincide with the name of any other gene.
    // impactful_mutations : is a comma-separated list of SNVs that need to be labelled with the value indicated in new_impact, format: chr5_1294991_C/T, with pyrimidine based definition
    // neutral_impact      : (optional, default; synonymous)
    // new_impact          : (optional, default: missense) is the impact that the mutations listed in impactful_mutations will receive.
    """
    panel_custom_processing.py \\
                    ${panel_annotated} \\
                    ${custom_regions} \\
                    ${panel_annotated.getBaseName()}.custom.tsv \\
                    ${simple} ;
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${panel_annotated.getBaseName()}.custom.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
