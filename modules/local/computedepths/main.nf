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
    path "versions.yml"                 , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def restrict_to_region = task.ext.restrict_panel ? "-b ${custombed}" : ""


    // this variable is used for subsetting the output depths table to only
    // positions with a mean depth above a given value
    // if the provided value is 0 this is not used
    def minimum_depth = task.ext.minimum_depth ? "| awk 'NR == 1 {print; next}  {sum = 0; for (i=3; i<=NF; i++) sum += \$i; mean = sum / (NF - 2); if (mean >= ${task.ext.minimum_depth} ) print }'": ""

    // this variable is used for setting a uniform depth across sites
    def uniform_depth = task.ext.uniform_depth ? "| awk 'NR == 1 {print; next} {for (i=3; i<=NF; i++) \$i = 100; print}'" : ""

    """
    ls -1 *.bam > bam_files_list.txt;
    samtools \\
        depth \\
        ${args} \\
        ${restrict_to_region} \\
        -@ $task.cpus \\
        -f bam_files_list.txt \\
        | tail -c +2 \\
        ${minimum_depth} \\
        ${uniform_depth} \\
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
