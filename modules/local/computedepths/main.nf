process COMPUTEDEPTHS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam)
    path (custombed)

    output:
    tuple val(meta), path("*.tsv.gz")   , emit: depths
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def restrict_to_region = task.ext.restrict_panel ? "-b ${custombed}" : ""
    """
    ls -1 *.bam > bam_files_list.txt;
    samtools \\
        depth \\
        ${args} \\
        ${restrict_to_region} \\
        -@ $task.cpus \\
        -f bam_files_list.txt \\
        | tail -c +2 \\
        | gzip -c > ${prefix}.depths.tsv.gz;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.depths.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}



// # Generate depth file using samtools depth
// samtools depth -H -@ $task.cpus -f bam_files_list.txt > ${prefix}.depths.tsv

// # Filter positions with mean depth >= 5 and compress the output
// awk '{
//     sum = 0;
//     for (i=3; i<=NF; i++) sum += $i;
//     mean = sum / (NF - 2);
//     if (mean >= 5) print
// }' ${prefix}.depths.tsv | gzip -c > ${prefix}.filtered.depths.tsv.gz

// # Optionally, remove the intermediate depth file if it's no longer needed
// rm ${prefix}.depths.tsv
