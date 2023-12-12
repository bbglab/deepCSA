process ONCODRIVEFML {
    // tag "$meta.id"
    tag "test"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/oncodrivefml:mutabilities'

    input:
    // tuple val(meta), path(mutations), path(mutabilities), path(bed_file)
    path(mutations)
    path(mutabilities)
    path(mutabilities_ind)
    path(bed_file)


    output:
    path("**.tsv.gz")  , emit: tsv
    path("**.png")     , emit: png
    path("**.html")    , emit: html
    // tuple val(meta), path("test/*") , emit: tsv
    // tuple val(meta), path("${prefix}/*.png")    , emit: plots
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "" // "-s ${params.seed}"
    def prefix = task.ext.prefix ?: "test" // "${meta.id}"
    def cadd_scores = task.ext.cadd_scores ?: '/workspace/datasets/CADD/v1.6/hg38/whole_genome_SNVs.tsv.gz'
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    cat > oncodrivefml_v2.mutability.conf << EOF
    [genome]
    build = 'hg38'

    [signature]
    method = 'none'


    [mutability]
    adjusting = True
    file = '$mutabilities'
    format = 'tabix'
    chr = 0
    chr_prefix = ""
    pos = 1
    ref = 2
    alt = 3
    mutab = 4


    [score]
    file = $cadd_scores
    format = 'tabix'
    chr = 0
    chr_prefix = ""
    pos = 1
    ref = 2
    alt = 3
    score = 5


    [statistic]
    method = 'amean'
    discard_mnp = False

    sampling = 100000
    sampling_max = 1000000
    sampling_chunk = 100
    sampling_min_obs = 10

        [[indels]]
            include = True
            method = 'max'
            max_consecutive = 7

    [settings]
    cores = $task.cpus
    EOF

    oncodrivefml -i $mutations \\
                    -e $bed_file \\
                    -o $prefix \\
                    $args \\
                    -c oncodrivefml_v2.mutability.conf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrivefml: \$(echo \$(oncodrivefml --version | rev | cut -d ' ' -f1 | rev))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "test" // "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrivefml: \$(echo \$(oncodrivefml --version | rev | cut -d ' ' -f1 | rev))
    END_VERSIONS
    """
}
