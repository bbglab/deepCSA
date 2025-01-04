process QUERY_BIOMART {

    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/dnds:latest'
    // we are probably missing python in this container

    input:
    tuple val(meta) , path(panel)
    tuple val(meta2), path(bed_file)

    output:
    tuple val(meta), path("custom_filtered_biomart.tsv"), emit: filtered_biomart
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f 5 ${panel} | sort -u | tr -s '\\n' ',' > genes_list.txt

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

    # Biomart Query
    biomart_url='http://jan2024.archive.ensembl.org/biomart/martservice'
    biomart_cds_query=\$(cat biomartQuery.txt)
    biomart_cds_query_encoded=\$(python -c "from urllib.parse import quote_plus; query = '''${biomart_cds_query}'''; print(quote_plus(query.replace('\\n', '')))")
    curl -L -s "${biomart_url}?query=\${biomart_cds_query_encoded}" | \\
        tail -n +2 | \\
        awk -F'\\t' '(\$5!=""){print(\$0)}' > biomart_output.tsv

    filter_biomart_query.R --bedfile ${bed_file} \\
                            --biomartoutput biomart_output.tsv \\
                            --outputfile custom_filtered_biomart.tsv


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

// "--referencetranscripts"
// default="/workspace/projects/prominent/analysis/dNdScv/data/reference_files/RefCDS_human_latest_intogen.rda",
// --covariates
//     "/workspace/projects/prominent/analysis/dNdScv/data/reference_files/covariates_hg19_hg38_epigenome_pcawg.rda",
//             help="Human GRCh38 covariates file [default= %default]", metavar="character"),
// --genelist"), type="character",
//             default=NULL,
//             help="Gene list file [default= %default]", metavar="character"),
// --genedepth"), type="character",
//             default=NULL,
//             help="Gene depth file (2 columns: GENE\tAVG_DEPTH) [default= %default]", metavar="character"),
// --snvsonly"), type="logical",
//             default=FALSE,
//             help="Only use SNVs for the analysis [default= %default]", metavar="logical")
    // <Filter name = "ensembl_gene_id" value = "ENSG00000141510,ENSG00000141511"/>
    // <Filter name = "external_gene_name" value = "`cat genes_list.txt`"/>
