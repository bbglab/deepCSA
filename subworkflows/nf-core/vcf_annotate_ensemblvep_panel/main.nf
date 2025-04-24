// taken from deepUMIcaller
//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/veppanel/main'
// include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'
include { SPLIT_TSV_BY_CHROM }    from '../../../modules/local/process_annotation/splittsvbychrom/main'
include { MERGE_VEP }      from '../../../modules/local/process_annotation/mergevep/main'



workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    vcf               // channel: [ val(meta), vcf ]
    fasta             //   value: fasta to use (optionnal)
    vep_genome        //   value: genome to use
    vep_species       //   value: species to use
    vep_cache_version //   value: cache version to use
    vep_cache         //    path: /path/to/vep/cache (optionnal)
    vep_extra_files   // channel: [ file1, file2...] (optionnal)

    main:

    SPLIT_TSV_BY_CHROM(tsv)
    // Run VEP on each chromosome in parallel
    ENSEMBLVEP_VEP(
        SPLIT_TSV_BY_CHROM.out.tsv_chunks.flatten(),
        vep_genome,
        vep_species,
        vep_cache_version,
        vep_cache,
        fasta,
        vep_extra_files
    )
    // Group the results by meta before merging
    vcf_annotated = ENSEMBLVEP_VEP.out.vcf.groupTuple()
    tab_annotated = ENSEMBLVEP_VEP.out.tab.groupTuple()
    json_annotated = ENSEMBLVEP_VEP.out.json.groupTuple()
    reports = ENSEMBLVEP_VEP.out.report.collect()

    // Merge the annotated tab files
    MERGE_VEP(tab_annotated)

    emit:
    vcf      = vcf_annotated      // channel: [ val(meta), [vcf1.gz, vcf2.gz, ...] ]
    tab      = MERGE_VEP.out.tab  // channel: [ val(meta), merged.tab.gz ]
    json     = json_annotated     // channel: [ val(meta), [json1.gz, json2.gz, ...] ]
    reports  = reports            // channel: [ *.html ]
}
