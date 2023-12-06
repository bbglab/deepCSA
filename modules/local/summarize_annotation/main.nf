process SUMMARIZE_ANNOTATION {
    // tag "$meta.id"
    tag "test"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    // tuple val(meta), path(mutations), path(mutabilities), path(bed_file)
    path(tab_files)

    output:
    path("*.summary.tab.gz")  , emit: tab
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "" // "-s ${params.seed}"
    def prefix = task.ext.prefix ?: "test" // "${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    zcat <\$(ls *.tab.gz | head -1) | grep '#' | grep -v '##' > header.tsv;
    for file in *.tab.gz; do
            zgrep -v '#' \$file >> all_samples.vep.tab.tmp;
            echo \$file;
    done;
    cat header.tsv <(sort -u all_samples.vep.tab.tmp | grep -v '#') > all_samples.vep.tab ;

    postprocessing_annotation.py all_samples.vep.tab all_samples.vep.summary.tab False;
    gzip all_samples.vep.summary.tab;

    rm all_samples.vep.tab.tmp;
    rm all_samples.vep.tab;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "test" // "${meta.id}"
    """
    touch all_samples.vep.summary.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

    // oncodrivefml -i $mutations \\
    //                 -e $bed_file \\
    //                 -o $prefix \\
    //                 $args \\
    //                 -c oncodrivefml_v2.mutability.conf
