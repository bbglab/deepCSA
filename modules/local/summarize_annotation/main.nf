process SUMMARIZE_ANNOTATION {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta), path(tab_files)

    output:
    tuple val(meta), path("*.summary.tab.gz")  , emit: tab
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assembly = task.ext.assembly ?: "hg38"

    // TODO reimplement python script with click
    """
    zcat <\$(ls *.tab.gz | head -1) | grep '#' | grep -v '##' > header.tsv;
    for file in *.tab.gz; do
            zgrep -v '#' \$file >> ${prefix}.vep.tab.tmp;
            echo \$file;
    done;
    cat header.tsv <(sort -u ${prefix}.vep.tab.tmp | grep -v '#') > ${prefix}.vep.tab ;
    rm ${prefix}.vep.tab.tmp;

    postprocessing_annotation.py ${prefix}.vep.tab ${prefix}.vep.summary.tab ${assembly} False;
    gzip ${prefix}.vep.summary.tab;

    rm ${prefix}.vep.tab;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
