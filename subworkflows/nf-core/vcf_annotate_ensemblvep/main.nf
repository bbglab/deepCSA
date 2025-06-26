//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'

// TODO update this so that it only receives the vcf as input
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

    ENSEMBLVEP_VEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)

    // Gather versions of all tools used

    emit:
    // vcf_tbi  = ch_vcf_tbi                  // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    vcf      = ENSEMBLVEP_VEP.out.vcf      // channel: [ val(meta), vcf ]
    json     = ENSEMBLVEP_VEP.out.json     // channel: [ val(meta), json ]
    tab      = ENSEMBLVEP_VEP.out.tab      // channel: [ val(meta), tab ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ *.html ]
}
