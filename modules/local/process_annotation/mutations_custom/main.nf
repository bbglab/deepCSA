process CUSTOM_MUTATION_PROCESSING {
    tag "$meta.id"

    label 'cpu_low'
    label 'process_high_memory'
    label 'time_low'


    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta)  , path(tab_files)
    path(hotspots_annotation_file)

    output:
    tuple val(meta), path("*.summary.tab.gz")   , emit: tab
    tuple val(meta), path("*.vep.tab.gz")       , emit: tab_all
    path "versions.yml"                         , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assembly = task.ext.assembly ?: "hg38"
    def hotspots = task.ext.hotspots_annotation ? "${hotspots_annotation_file}" : ""
    // TODO reimplement python script with click
    """
    zcat <\$(ls *.tab.gz | head -1) | grep '#' | grep -v '##' > header.tsv;
    for file in *.tab.gz; do
            zgrep -v '#' \$file >> ${prefix}.vep.tab.tmp;
            echo \$file;
    done;
    cat header.tsv <(sort -u ${prefix}.vep.tab.tmp | grep -v '#' | sed 's/^chr//g' | awk -F ':' 'length(\$1) <= 2' | awk '{ printf "chr"\$0"\\n" }') > ${prefix}.vep.tab ;
    rm ${prefix}.vep.tab.tmp;

    postprocessing_annotation.py ${prefix}.vep.tab \\
                                    ${prefix}.vep.summary.tab \\
                                    ${assembly} \\
                                    False \\
                                    ${hotspots} ;
    gzip ${prefix}.vep.summary.tab;
    gzip ${prefix}.vep.tab;

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
