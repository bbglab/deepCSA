process ONCODRIVEFML {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/oncodrivefml:latest'

    input:
    tuple val(meta) , path(mutations), path(mutabilities), path(mutabilities_ind)
    tuple val(meta2), path (bed_file)
    tuple path(cadd_scores_file), path(cadd_scores_index)
    val(mode)

    output:
    tuple val(meta), path("**.tsv.gz")  , emit: tsv
    tuple val(meta), path("**.png")     , emit: png
    tuple val(meta), path("**.html")    , emit: html
    tuple val(meta), path("${meta.id}.${mode}/")         , emit: folder
    path "versions.yml"                 , topic: versions


    script:
    def args = task.ext.args ?: "" // "-s ${params.seed}"
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    // TODO: See if we can provide the entire json as an input parameter
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
    chr_prefix = "chr"
    pos = 1
    ref = 2
    alt = 3
    mutab = 4


    [depth]
    adjusting = False
    file = '${mutabilities}'
    format = 'tabix'
    chr = 0
    chr_prefix = "chr"
    pos = 1
    depth = 2


    [score]
    file = ${cadd_scores_file}
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
            gene_exomic_frameshift_ratio = True
            stops_function = 'random_choice'

    [settings]
    cores = ${task.cpus}
    seed = 123
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
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrivefml: \$(echo \$(oncodrivefml --version | rev | cut -d ' ' -f1 | rev))
    END_VERSIONS
    """
}
