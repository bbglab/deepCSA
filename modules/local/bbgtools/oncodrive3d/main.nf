process ONCODRIVE3D {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    // TODO pending to push the container somewhere and be able to retrieve it
    container 'docker.io/ferriolcalvet/oncodrive3d:latest'

    input:
    tuple val(meta), path(mutations), path(mutabilities), path(mutabilities_ind)
    path (datasets)


    output:
    tuple val(meta), path("**genes.csv")                 , emit: csv_genes
    tuple val(meta), path("**pos.csv")                   , emit: csv_pos
    tuple val(meta), path("**mutations.processed.tsv")   , emit: mut_processed
    tuple val(meta), path("**miss_prob.processed.json")  , emit: prob_processed
    tuple val(meta), path("**seq_df.processed.tsv")      , emit: seq_processed
    tuple val(meta), path("**.log")                      , emit: log
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mane = false
    def vep_raw = false
    """
    cat > oncodrive3d.mutability.conf << EOF
    {
        "file" : "${mutabilities}",
        "format" : "tabix",
        "chr" : 0,
        "chr_prefix" : "chr",
        "pos" : 1,
        "ref" : 2,
        "alt" : 3,
        "mutab" : 4
    }
    EOF

    oncodrive3D run -i ${mutations} \\
                    -m oncodrive3d.mutability.conf \\
                    -d ${datasets} \\
                    -C ${prefix} \\
                    -o ${prefix} \\
                    ${args} \\
                    -c ${task.cpus} \\
                    ${vep_raw ? '--o3d_transcripts --use_input_symbols' : ''} \\
                    ${mane ? '--mane' : ''} \\
                    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: 1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: 2.0
    END_VERSIONS
    """
}

