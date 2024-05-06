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
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tool = task.ext.tool ?: 'oncodrivefml'

    // find a better solution for doing this,
    // probably in python
    // so that both genes overlapping can be conserved and then we can also do groups of genes
    // step 1 make it robust
    // step 2 allow groups
    """
    if [ "$tool" == "oncodrivefml" ]; then
        cat <(printf "CHROMOSOME\\tSTART\\tEND\\tELEMENT\\tSEGMENT\\n") <(cut -f1,2,6,7 ${panel_tsv} | tail -n +2 | uniq | awk -F'\\t' '{print \$1, \$2, \$2, \$3}' OFS='\\t' | bedtools merge -c 4 -o distinct | sed 's/^chr//g' | awk -F'\\t' '{print \$1, \$2, \$3, \$4, 1}' OFS='\\t' ) > ${prefix}.annotated.bed
    elif [ "$tool" == "oncodriveclustl" ]; then
        cat <(printf "CHROMOSOME\\tSTART\\tEND\\tELEMENT_ID\\tSYMBOL\\n") <(cut -f1,2,6,7 ${panel_tsv} | tail -n +2 | uniq | awk -F'\\t' '{print \$1, \$2, \$2, \$3}' OFS='\\t' | bedtools merge -c 4 -o distinct | sed 's/^chr//g' | awk -F'\\t' '{print \$1, \$2, \$3, \$4, \$4}' OFS='\\t' ) > ${prefix}.annotated.bed
    elif [ "$tool" == "readsxposition" ]; then
        cat <(printf "CHROMOSOME\\tSTART\\tEND\\tGENE\\tSYMBOL\\n") <(cut -f1,2,6,7 ${panel_tsv} | tail -n +2 | uniq | awk -F'\\t' '{print \$1, \$2, \$2, \$3}' OFS='\\t' | bedtools merge -c 4 -o distinct | awk -F'\\t' '{print \$1, \$2, \$3, \$4, \$4}' OFS='\\t' ) > ${prefix}.annotated.bed
    fi

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
