process ENSEMBLVEP_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    conda params.vep_cache_version == 108 ? 'bioconda::ensembl-vep=108.2' : 
            params.vep_cache_version == 102 ? 'bioconda::ensembl-vep=102.0' :  
            params.vep_cache_version == 111 ? 'bioconda::ensembl-vep=111.0' :  
            'bioconda::ensembl-vep=111.0' 

    container params.vep_cache_version == 108 ? "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/ensembl-vep:108.2--pl5321h4a94de4_0' : 'biocontainers/ensembl-vep:108.2--pl5321h4a94de4_0' }" : 
                params.vep_cache_version == 102 ? "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/ensembl-vep:102.0--pl526hecda079_0' : 'biocontainers/ensembl-vep:102.0--pl526hecda079_0' }" : 
                params.vep_cache_version == 111 ? "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/ensembl-vep:111.0--pl5321h2a3209d_0' : 'biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0' }" :
                "biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0"

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)

    output:
    tuple val(meta), path("vep_cache"), emit: cache
    path "versions.yml"               , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vep_install \\
        --CACHEDIR vep_cache \\
        --SPECIES $species \\
        --ASSEMBLY $assembly \\
        --CACHE_VERSION $cache_version \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir vep_cache

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}