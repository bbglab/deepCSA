#!/usr/bin/env python

# if the maf is from deepCSA, use this one, otherwise use the one below
# usage: python deepcsa_maf2vcf.py --mutations-file all_samples.somatic.mutations.tsv --output-dir ./test/ --maf-from-deepcsa

# if the maf file is not from deepCSA, use this below
# usage: python deepcsa_maf2vcf.py --mutations-file all_samples.somatic.mutations.tsv --output-dir ./test/

import click
import pandas as pd


def build_vcf_like_dataframe(mutations_dataframe, samplee):
    """
    Build a VCF-like dataframe from the mutations dataframe.
    input needs to have:
        ['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'ALT_DEPTH']
    output needs to have:
        ['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
    """
    for col in ['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'ALT_DEPTH']:
        if col not in mutations_dataframe.columns:
            raise ValueError(f"Column {col} is missing from the mutations dataframe.")
    
    # fill FILTER and INFO columns with default values
    if "FILTER" not in mutations_dataframe.columns:
        print("WARNING: FILTER column is missing from the mutations dataframe. Setting it to 'PASS' for all mutations")
        mutations_dataframe["FILTER"] = "PASS"
    if "INFO" not in mutations_dataframe.columns:
        print(f"WARNING: INFO column is missing from the mutations dataframe. Setting it to 'SAMPLE={samplee};'")
        mutations_dataframe["INFO"] = f"SAMPLE={samplee};"

    # Create a new dataframe with the required columns
    vcf_like_df = mutations_dataframe[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO', 'DEPTH', 'ALT_DEPTH']].copy()
    vcf_like_df["FORMAT"] = "GT:DP:VD:AD:AF:RD:ALD:CDP:CAD:NDP:CDPAM:CADAM:NDPAM"
    vcf_like_df["SAMPLE"] = vcf_like_df[['DEPTH', 'ALT_DEPTH']].apply(
        lambda x: "{GT}:{DP}:{VD}:{AD}:{AF}:{RD}:{ALD}:{CDP}:{CAD}:{NDP}:{CDPAM}:{CADAM}:{NDPAM}".format(
            GT="0/1",
            DP=x['DEPTH'],
            VD=x['ALT_DEPTH'],
            AD=f"{x['DEPTH'] - x['ALT_DEPTH']},{x['ALT_DEPTH']}",
            AF=round(x['ALT_DEPTH'] / x['DEPTH'] , 5),
            RD=f"{(x['DEPTH'] - x['ALT_DEPTH'])//2},{(x['DEPTH'] - x['ALT_DEPTH'])//2 if (x['DEPTH'] - x['ALT_DEPTH']) % 2 == 0 else (x['DEPTH'] - x['ALT_DEPTH'])//2 + 1}",
            ALD=f"{x['ALT_DEPTH']//2},{x['ALT_DEPTH']//2 if x['ALT_DEPTH'] % 2 == 0 else x['ALT_DEPTH']//2 + 1}",
            CDP=x['DEPTH'],
            CAD=f"{x['DEPTH'] - x['ALT_DEPTH']},{x['ALT_DEPTH']}",
            NDP="0",
            CDPAM=x['DEPTH'],
            CADAM=f"{x['DEPTH'] - x['ALT_DEPTH']},{x['ALT_DEPTH']}",
            NDPAM="0" 
        ),
        axis=1
    )
    return vcf_like_df

filters_to_remove = ["not_in_exons", "not_covered"]
def remove_deepcsa_filters(old_filt, filters_to_removee):
    """
    Remove deepCSA filters from the FILTER field of the VCF file.
    """    
    filter_result = sorted([ x for x in old_filt.split(";") if x not in filters_to_removee ])
    return ";".join(filter_result) if filter_result != [] else "PASS"


vardict_vcf_header = '''##fileformat=VCFv4.2
##source=artificially_generated_vcf_VarDict_v1.8.2-like
##FILTER=<ID=AMPBIAS,Description="Indicate the variant has amplicon bias.">
##FILTER=<ID=AM_no_pileup_support,Description="Variant not supported when inspecting the BAM with mpileup. deepUMIcaller.">
##FILTER=<ID=AM_not_searched_COMPLEX,Description="Complex variant. Recounting of alt depth from mpileup output not done for this type of variants. deepUMIcaller.">
##FILTER=<ID=AM_not_searched_SV,Description="Structural variant. Recounting of alt depth from mpileup output not done for this type of variants. deepUMIcaller.">
##FILTER=<ID=Bias,Description="Strand Bias">
##FILTER=<ID=Cluster0bp,Description="Two variants are within 0 bp">
##FILTER=<ID=InGap,Description="The variant is in the deletion gap, thus likely false positive">
##FILTER=<ID=InIns,Description="The variant is adjacent to an insertion variant">
##FILTER=<ID=LongMSI,Description="The somatic variant is flanked by long A/T (>=14)">
##FILTER=<ID=MSI12,Description="Variant in MSI region with 12 non-monomer MSI or 13 monomer MSI">
##FILTER=<ID=NM20,Description="Mean mismatches in reads >= 20, thus likely false positive">
##FILTER=<ID=Q10,Description="Mean Mapping Quality Below 10">
##FILTER=<ID=SN1.5,Description="Signal to Noise Less than 1.5">
##FILTER=<ID=d3,Description="Total Depth < 3">
##FILTER=<ID=f0.0,Description="Allele frequency < 0.0">
##FILTER=<ID=low_complex_repetitive,Description="Variant located in a low complexity or repetitive genomic region. deepUMIcaller.">
##FILTER=<ID=low_mappability,Description="Variant located in a genomic region with low mappability. deepUMIcaller.">
##FILTER=<ID=n_rich,Description="Variant located in a position with more Ns than expected. Threshold in this sample: 0.491005. deepUMIcaller.">
##FILTER=<ID=no_pileup_support,Description="Variant not supported when inspecting the BAM with mpileup. deepUMIcaller.">
##FILTER=<ID=not_searched_COMPLEX,Description="Complex variant. Recounting of alt depth from mpileup output not done for this type of variants. deepUMIcaller.">
##FILTER=<ID=not_searched_SV,Description="Structural variant. Recounting of alt depth from mpileup output not done for this type of variants. deepUMIcaller.">
##FILTER=<ID=p8,Description="Mean Position in Reads Less than 8">
##FILTER=<ID=pSTD,Description="Position in Reads has STD of 0">
##FILTER=<ID=q22.5,Description="Mean Base Quality Below 22.5">
##FILTER=<ID=v2,Description="Var Depth < 2">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">
##FORMAT=<ID=CAD,Number=2,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. ONLY ONE ALT allele. Recomputed using mpileup output, this might not sum to CDP. deepUMIcaller.">
##FORMAT=<ID=CADAM,Number=2,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. ONLY ONE ALT allele. Recomputed using mpileup output, this might not sum to CDP. deepUMIcaller.">
##FORMAT=<ID=CDP,Number=1,Type=Integer,Description="Total Depth recomputed using mpileup output. deepUMIcaller.">
##FORMAT=<ID=CDPAM,Number=1,Type=Integer,Description="Total Depth recomputed using mpileup output. deepUMIcaller.">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=NDP,Number=1,Type=Integer,Description="Total number of Ns computed using mpileup output. deepUMIcaller.">
##FORMAT=<ID=NDPAM,Number=1,Type=Integer,Description="Total number of Ns computed using mpileup output. deepUMIcaller.">
##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##INFO=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AMPFLAG,Number=1,Type=Integer,Description="Top variant in amplicons don't match">
##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand Bias Info">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=DUPRATE,Number=1,Type=Float,Description="Duplication rate in fraction">
##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">
##INFO=<ID=GDAMP,Number=1,Type=Integer,Description="No. of amplicons supporting variant">
##INFO=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">
##INFO=<ID=HICNT,Number=1,Type=Integer,Description="High quality variant reads">
##INFO=<ID=HICOV,Number=1,Type=Integer,Description="High quality total reads">
##INFO=<ID=LSEQ,Number=1,Type=String,Description="5' flanking seq">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">
##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSatellite. > 1 indicates MSI">
##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MicroSatellite unit length in bp">
##INFO=<ID=NCAMP,Number=1,Type=Integer,Description="No. of amplicons don't work">
##INFO=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">
##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Odds ratio">
##INFO=<ID=PMEAN,Number=1,Type=Float,Description="The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call">
##INFO=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">
##INFO=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">
##INFO=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">
##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Reference depth by strand">
##INFO=<ID=RSEQ,Number=1,Type=String,Description="3' flanking seq">
##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">
##INFO=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">
##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">
##INFO=<ID=SN,Number=1,Type=Float,Description="Signal to noise">
##INFO=<ID=SPANPAIR,Number=1,Type=Integer,Description="No. of pairs supporting SV">
##INFO=<ID=SPLITREAD,Number=1,Type=Integer,Description="No. of split reads supporting SV">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="The length of SV in bp">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type: INV DUP DEL INS FUS">
##INFO=<ID=TLAMP,Number=1,Type=Integer,Description="Total of amplicons covering variant">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">
##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Variant depth by strand">
##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
'''


@click.command()
@click.option('--mutations-file', required=True, type=click.Path(exists=True), help="Path to the mutations file (TSV format).")
@click.option('--output-dir', required=True, type=click.Path(), help="Directory to save the output VCF files.")
@click.option('--maf-from-deepcsa', is_flag=True, default=False, help="Flag to indicate if the MAF file is from deepCSA.")
def main(mutations_file, output_dir, maf_from_deepcsa):
    """
    Convert a mutations file to one or multiple VCF-formatted files.
    """

    mutations = pd.read_table(mutations_file)

    for sample in mutations["SAMPLE_ID"].unique():
        
        # filter the dataframe for the current sample
        sample_mutations = mutations[mutations["SAMPLE_ID"] == sample]

        if maf_from_deepcsa:
            vcf_info_sample = sample_mutations[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']].copy()
        else:
            vcf_info_sample = build_vcf_like_dataframe(sample_mutations, sample)
            # mandatory columns should be: [['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'DEPTH', 'ALT_DEPTH']]
            # Ns can be assumed 0 and AM can be assumed to be the same as duplex

        # clean FILTER field of all deepCSA annotations
        if len(filters_to_remove) > 0:
            vcf_info_sample["FILTER"] = vcf_info_sample["FILTER"].apply(remove_deepcsa_filters, filters_to_removee = filters_to_remove)
        
        # add other necessary columns
        vcf_info_sample["ID"] = '.'
        vcf_info_sample["QUAL"] = 100

        # rename sample column
        vcf_info_sample.rename(columns={"SAMPLE": sample}, inplace=True)
        
        cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]
        
        sample_file = f"{output_dir}/{sample}.vcf"
        
        with open(sample_file, "w") as f:
            f.write(vardict_vcf_header)
            f.write("#" + "\t".join(cols) + "\n")
            
            # write the data
            vcf_info_sample.to_csv(f, sep="\t", columns=cols, index=False, header=False)
            
        print(f"VCF file for {sample} : {sample_file}")

if __name__ == "__main__":
    main()