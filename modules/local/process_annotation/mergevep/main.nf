process MERGE_VEP {
    label 'process_medium'

    conda "bioconda::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/pigz:2.3.4' : 'quay.io/biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(vep_files)

    output:
    tuple val(meta), path("*.merged.tab.gz"), emit: tab
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Output filename
    out_file="${prefix}.merged.tab"

    # Sort input files by chromosome (e.g. chr1, chr2, ..., chrX, chrY)
    sorted_files=\$(printf "%s\\n" ${vep_files} | sort -V)

    # Initialize header_written flag
    header_written=false

    # Merge files
    for f in \$sorted_files; do
        echo "Processing \$f..."
        if [ "\$header_written" = false ]; then
            # Use pigz for parallel decompression and awk for header and data processing
            pigz -dc -p ${task.cpus} "\$f" | awk '{if(\$0 ~ /^##/) print; else if(\$0 ~ /^#/) {print; next} else print}' > "\$out_file"
            header_written=true
        else
            # Only append data rows (skip header lines)
            pigz -dc -p ${task.cpus} "\$f" | awk 'NR>1 {print}' >> "\$out_file"
        fi
    done

    # Compress final merged output using pigz instead of bgzip
    pigz -p ${task.cpus} "\$out_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.merged.tab.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tools: gzip \$(gzip --version | head -n1)
    END_VERSIONS
    """
}