process ONCODRIVECLUSTL {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/oncodriveclustl:latest'

    input:
    tuple val(meta) , path (mutations), path(mutabilities), path(mutabilities_ind)
    tuple val(meta2), path (bed_file)

    output:
    tuple val(meta), path("**.tsv")                 , emit: tsv
    tuple val(meta), path("**.txt")                 , emit: txt
    tuple val(meta), path("**.png")  , optional:true, emit: png
    tuple val(meta), path("**.log")                 , emit: log
    path "versions.yml"                             , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def assembly = task.ext.assembly ?: "hg38"
    """
    cat > mutability_config.json << EOF
    {
        "file" : "$mutabilities",
        "format" : "tabix",
        "chr" : 0,
        "chr_prefix" : "chr",
        "pos" : 1,
        "ref" : 2,
        "alt" : 3,
        "mutab" : 4
    }
    EOF

    oncodriveclustl -i ${mutations} \\
                    -r ${bed_file} \\
                    -o ${prefix} \\
                    --cores ${task.cpus} \\
                    -g ${assembly} \\
                    ${args} \\
                    -mutab mutability_config.json

    # cd ~/projects/oncodriveclustl_dev/oncodriveclustl
    # oncodriveclustl -i example/hg38_test2/mutations.all_samples.tsv.gz -r example/hg38_test2/KidneyGenes.canonical_transcripts_CDS.expanded2.in_panel.bed6.oncodriveclustl.bed -o example/hg38_test2/results_mutability_mut_centered -g hg38 -sim mutation_centered -kmer 3 -mutab example/hg38_test2/mutability_config.json --clustplot


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodriveclustl: \$(echo \$(oncodriveclustl --version | rev | cut -d ' ' -f1 | rev))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodriveclustl: \$(echo \$(oncodriveclustl --version | rev | cut -d ' ' -f1 | rev))
    END_VERSIONS
    """
}
