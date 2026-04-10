process QUERY_BIOMART {

    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    label 'deepcsa_core'

    input:
    tuple val(meta) , path(panel)
    tuple val(meta2), path(bedfile)

    output:
    tuple val(meta), path("custom_filtered_biomart.tsv"), emit: filtered_biomart
    tuple val(meta), path("splice_sites.tsv")           , emit: splice_sites
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cut -f 6 ${panel} | tail -n +2 | sort -u | awk '\$1!="-"' | tr -s '\\n' ',' > genes_list.txt

    cat > biomartQuery.txt << EOF
    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
        <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
            <Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
            <Filter name = "biotype" value = "protein_coding"/>
            <Filter name = "external_gene_name" value = "\$(cat genes_list.txt)"/>
            <Filter name = "mane_select" excluded = "0"/>
            <Attribute name = "ensembl_gene_id" />
            <Attribute name = "external_gene_name" />
            <Attribute name = "ensembl_peptide_id" />
            <Attribute name = "chromosome_name" />
            <Attribute name = "genomic_coding_start" />
            <Attribute name = "genomic_coding_end" />
            <Attribute name = "cds_start" />
            <Attribute name = "cds_end" />
            <Attribute name = "cds_length" />
            <Attribute name = "strand" />
            <Attribute name = "ensembl_transcript_id" />
            <Attribute name = "exon_chrom_start" />
            <Attribute name = "exon_chrom_end" />
        </Dataset>
    </Query>
    EOF

    wget -O result.txt \
        --post-file=query.xml \
        'https://jan2024.archive.ensembl.org/biomart/martservice'

    dNdScv_panel_prep.py --bed ${bedfile} \\
                        --genes result.txt \\
                        --output custom_filtered_biomart.tsv \\
                        --verbose \
                        -s splice_sites.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ensembl BioMart: v111
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.out.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ensembl BioMart: v111
    END_VERSIONS
    """
}