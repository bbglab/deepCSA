// taken from deepUMIcaller

process ENSEMBLVEP_VEP {
    tag "$meta.id"
    label 'process_high'
    label 'cpu_veryhigh'

    conda params.vep_cache_version == 108 ? 'bioconda::ensembl-vep=108.2' : 
            params.vep_cache_version == 102 ? 'bioconda::ensembl-vep=102.0' :  
            params.vep_cache_version == 111 ? 'bioconda::ensembl-vep=111.0' :  
            'bioconda::ensembl-vep=111.0' 

    container params.vep_cache_version == 108 ? "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/ensembl-vep:108.2--pl5321h4a94de4_0' : 'biocontainers/ensembl-vep:108.2--pl5321h4a94de4_0' }" : 
                params.vep_cache_version == 102 ? "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/ensembl-vep:102.0--pl526hecda079_0' : 'biocontainers/ensembl-vep:102.0--pl526hecda079_0' }" : 
                params.vep_cache_version == 111 ? "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'https://depot.galaxyproject.org/singularity/ensembl-vep:111.0--pl5321h2a3209d_0' : 'biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0' }" :
                "biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0"


    input:
    tuple val(meta), path(vcf)
    val   genome
    val   species
    val   cache_version
    path  cache
    path  fasta
    path  extra_files

    output:
    tuple val(meta), path("*.vcf.gz")  , optional:true, emit: vcf
    tuple val(meta), path("*.tab.gz")  , optional:true, emit: tab
    tuple val(meta), path("*.json.gz") , optional:true, emit: json
    path "*.summary.html"              , optional:true, emit: report
    path "versions.yml"                               , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""

    """

    # Convert input TSV to VEP format, to make vep --fork more efficient
    awk 'BEGIN { OFS="\t" }
    {
     split(\$4, a, "/");
     print \$1, \$2, ".", a[1], a[2];
    }' ${vcf} > ${vcf}.vep

    # this is to ensure that we will be able to match the tab and vcf files afterwards
    # the structure of the ID is the following:
    vep \\
        -i ${vcf}.vep \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf.gz
    touch ${prefix}.ann.tab.gz
    touch ${prefix}.ann.json.gz
    touch ${prefix}.summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
