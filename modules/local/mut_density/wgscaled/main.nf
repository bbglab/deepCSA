process WG_SCALED_MUTATION_DENSITY {

    tag "${meta.id}"

    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/runningr:v1'

    input:
    tuple val(meta), path(mutations_file), path(depths_file)
    tuple val(meta2), path(consensus_bed_all)
    path (wgs_counts)

    output:
    tuple val(meta), path("*.tsv")  , emit: adjusted_mutrate
    path "versions.yml"             , topic: versions


    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    mutrate_genome_trinuc_corrected.R \\
        --samplename ${meta.id} \\
        --mutations ${mutations_file} \\
        --depths ${depths_file} \\
        --consensus_bed ${consensus_bed_all} \\
        --wgs_counts ${wgs_counts} \\
        --outputprefix ${prefix} \\
        --panel_version ${panel_version}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.3.1
        Rscript: 4.3.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}_mutrates_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.3.1
        Rscript: 4.3.1
    END_VERSIONS
    """
}