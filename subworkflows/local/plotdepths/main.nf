include { TABIX_BGZIPTABIX_QUERY    as QUERYDEPTHS     } from '../../../modules/nf-core/tabix/bgziptabixquery/main'

include { PLOT_DEPTHS               as DEPTHSSUMMARY    } from '../../../modules/local/plot/depths_summary/main'

include { CREATECUSTOMBEDFILE       as CUSTOMBEDFILE    } from '../../../modules/local/createpanels/custombedfile/main'


workflow PLOT_DEPTHS {

    take:
    depth
    bedfile
    panel

    main:
    // Define the chromosome region
    chromosome_region = params.chromosome_region_file
                            ? channel.fromPath(params.chromosome_region_file, checkIfExists: true).first()
                            : channel.fromPath(file("${projectDir}/assets/placeholder_no_file.tsv", checkIfExists: true)).first()

    // Intersect BED of all sites with BED of sample filtered sites
    QUERYDEPTHS(depth, bedfile)

    CUSTOMBEDFILE(panel)

    DEPTHSSUMMARY(QUERYDEPTHS.out.subset, CUSTOMBEDFILE.out.bed, chromosome_region)


    emit:
    plots                       = DEPTHSSUMMARY.out.plots
    average_depth_sample        = DEPTHSSUMMARY.out.average_per_sample
    average_depth_gene          = DEPTHSSUMMARY.out.average_per_gene
    average_depth_gene_sample   = DEPTHSSUMMARY.out.average_per_gene_sample
    average_depth_region        = DEPTHSSUMMARY.out.average_per_region
    average_depth_region_sample = DEPTHSSUMMARY.out.average_per_region_sample
    genes_per_region            = DEPTHSSUMMARY.out.genes_per_region

}
