process CREATECUSTOMBEDFILE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pybedtools=0.9.1--py38he0f268d_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/pybedtools:0.9.1--py38he0f268d_0' :
            'biocontainers/pybedtools:0.9.1--py38he0f268d_0' }"

    input:
    tuple val(meta) , path(panel_tsv)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tool = task.ext.tool ?: 'oncodrivefml'
    // find a better solution for doing this,
    // probably in python
    // so that both genes overlapping can be conserved and then we can also do groups of genes
    // step 1 make it robust
    // step 2 allow groups
    """
    sh createcustombed.sh ${panel_tsv} ${tool} > ${prefix}.annotated.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotated.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
