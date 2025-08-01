process ONCODRIVE3D_RUN {
    tag "$meta.id"
    label 'process_high'

    container 'docker.io/bbglab/oncodrive3d:1.0.5'


    input:
    tuple val(meta), path(mutations), path(mutabilities), path(mutabilities_ind)
    path (datasets)


    output:
    tuple val(meta), path("**genes.csv")                 , emit: csv_genes
    tuple val(meta), path("**pos.csv")                   , emit: csv_pos
    tuple val(meta), path("**mutations.processed.tsv")   , emit: mut_processed, optional: true
    tuple val(meta), path("**miss_prob.processed.json")  , emit: prob_processed, optional: true
    tuple val(meta), path("**seq_df.processed.tsv")      , emit: seq_processed, optional: true
    tuple val(meta), path("**run_*.log")                 , emit: log
    path "versions.yml"                                  , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    // To be set as true to prioritize MANE transcripts but penalize plotting
    // annotations (we should also change the datasets dir to the MANE one)
    def mane = task.ext.mane ? '--mane' : ''
    def vep_raw = task.ext.vep_raw ? '--o3d_transcripts --use_input_symbols' : ''
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

    oncodrive3D run -i $mutations \\
                    -m oncodrive3d.mutability.conf \\
                    -d $datasets \\
                    -C $prefix \\
                    -o $prefix \\
                    $args \\
                    -c ${task.cpus} \\
                    ${vep_raw} \\
                    ${mane}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: \$(oncodrive3d --version | rev | cut -d ' ' -f 1 | rev )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrive3D: \$(oncodrive3d --version | rev | cut -d ' ' -f 1 | rev )
    END_VERSIONS
    """
}

