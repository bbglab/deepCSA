process SITESFROMPOSITIONS {

    tag "${meta.id}"

    container "docker.io/bbglab/deepcsa-core:0.0.2-alpha"


    input:
    tuple val(meta), path(depths)

    output:
    tuple val(meta), path("*.sites4VEP.chunk*.tsv")  , emit: annotated_panel_reg
    path  "versions.yml"                             , topic: versions


    script:
    def assembly = task.ext.assembly ?: "hg38"
    def chunk_size = task.ext.chunk_size ?: 0

    // TODO
    // see if there is a better way to filter out chromosomes
    // that are not the canonical ones right now doing it by
    // filtering the size of the chromosome name to be smaller of equal to 2

    """
    cat <(printf "CHROM\\tPOS\\n") <( zcat ${depths} | cut -f1,2  | sed 's/^chr//g' | awk 'length(\$1) <= 2' ) > captured_positions.tsv;
    sites_table_from_positions.py \\
        --input-positions captured_positions.tsv \\
        --genome-assembly ${assembly} \\
        --output-file-with-sites captured_positions.sites4VEP.tmp.tsv;

    rm captured_positions.tsv

    awk '{print "chr"\$0}' captured_positions.sites4VEP.tmp.tsv > captured_positions.sites4VEP.full.tsv

    # Chunk the sites file if chunk_size is set
    if [ ${chunk_size} -gt 0 ]; then
        echo "[SITESFROMPOSITIONS] Chunking sites file with chunk_size=${chunk_size}"
        
        # Extract header
        head -n 1 captured_positions.sites4VEP.full.tsv > header.tmp
        
        # Split file into chunks (excluding header)
        tail -n +2 captured_positions.sites4VEP.full.tsv | split -l ${chunk_size} --additional-suffix=.tsv -d - captured_positions.sites4VEP.chunk
        
        # Add header to each chunk
        for chunk in captured_positions.sites4VEP.chunk*.tsv; do
            cat header.tmp "\$chunk" > "\${chunk}.tmp" && mv "\${chunk}.tmp" "\$chunk"
        done
        
        n_chunks=\$(ls captured_positions.sites4VEP.chunk*.tsv | wc -l)
        echo "[SITESFROMPOSITIONS] Created \${n_chunks} chunks"
        
        rm header.tmp captured_positions.sites4VEP.full.tsv
    else
        echo "[SITESFROMPOSITIONS] No chunking (chunk_size=0), processing as single file"
        mv captured_positions.sites4VEP.full.tsv captured_positions.sites4VEP.chunk1.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch captured_positions.sites4VEP.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
