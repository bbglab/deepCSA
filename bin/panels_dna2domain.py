#!/usr/bin/env python


import click
import pandas as pd



def parse_panel_n_domains(consensus_panel_rich_file, domains_file):
    """
    Parse the consensus panel rich file and domains file to extract relevant information.
    The domains file must have these four columns:
        'Ens_Transcr_ID', 'Begin', 'End', 'NAME'
    """
    all_protein_positions = pd.read_table(consensus_panel_rich_file)[['CHROM', 'POS', 'STRAND', 'GENE', 'Feature', 'Protein_position']]
    all_protein_positions = all_protein_positions[all_protein_positions["Protein_position"] != "-"].drop_duplicates()
    all_protein_positions["Protein_position"] = all_protein_positions["Protein_position"].astype(float)
    all_protein_positions_summary = all_protein_positions[~(all_protein_positions["Protein_position"].isna())
                                                            ][['GENE', 'Protein_position', "CHROM", "POS", "STRAND"]].reset_index(drop = True)

    panel_all_gene_transcript = all_protein_positions[["GENE", "Feature"]].drop_duplicates()

    pfam_domains = pd.read_table(domains_file, usecols= ['Ens_Transcr_ID', 'Begin', 'End', 'NAME']).drop_duplicates()
    pfam_domains_summary = pfam_domains.merge(panel_all_gene_transcript,
                                                left_on = 'Ens_Transcr_ID',
                                                right_on = 'Feature',
                                                how = 'inner'
                                                ).drop(columns = ['Feature'])
    pfam_domains_summary["NAME"] = pfam_domains_summary["GENE"] + '--' + pfam_domains_summary["NAME"]

    return all_protein_positions_summary, pfam_domains_summary




def generate_domains2dna_mapping(consensus_panel_rich_file, domains_file, output_file):
    """
    Main function to process DNA to protein mapping.
    """

    protein_positions_summary, domains_info = parse_panel_n_domains(consensus_panel_rich_file, domains_file)

    all_protein_positions_small_summary = protein_positions_summary.sort_values(by = ["CHROM", "POS"],
                                                                                    ascending = True).drop_duplicates(
                                                                                                            subset = ['GENE', 'Protein_position'],
                                                                                                            keep = 'first')
    all_protein_positions_small_summary.columns = ['GENE', 'PROT_POS', 'CHROM', 'START_fw', "STRAND"]
    all_protein_positions_big_summary = protein_positions_summary.sort_values(by = ["CHROM", "POS"],
                                                                                    ascending = True).drop_duplicates(
                                                                                                            subset = ['GENE', 'Protein_position'],
                                                                                                            keep = 'last')
    all_protein_positions_big_summary.columns = ['GENE', 'PROT_POS', 'CHROM', 'END_fw', "STRAND"]



    pfam_domains_genomic = domains_info.merge(all_protein_positions_small_summary,
                                                                left_on = ['GENE', 'Begin'],
                                                                right_on = ["GENE", 'PROT_POS'],
                                                                how = 'left').merge(all_protein_positions_big_summary,
                                                                left_on = ['GENE', 'End', "CHROM", "STRAND"],
                                                                right_on = ["GENE", 'PROT_POS', "CHROM", "STRAND"],
                                                                how = 'left')
    pfam_domains_genomic_summary = pfam_domains_genomic[['CHROM', 'START_fw', 'END_fw', 'NAME', "STRAND"]].copy()
    pfam_domains_genomic_summary = pfam_domains_genomic_summary[
                                            (~(pfam_domains_genomic_summary["CHROM"].isna())) &
                                            (~(pfam_domains_genomic_summary["START_fw"].isna())) &
                                            (~(pfam_domains_genomic_summary["END_fw"].isna()))
                                        ].reset_index(drop = True)
    pfam_domains_genomic_summary["START_fw"] = pfam_domains_genomic_summary["START_fw"].astype(int)
    pfam_domains_genomic_summary["END_fw"] = pfam_domains_genomic_summary["END_fw"].astype(int)



    pfam_domains_genomic_summary[["START","END"]] =  pfam_domains_genomic_summary[["START_fw", "END_fw", "STRAND"]].apply(
                                                                    lambda x : (x["END_fw"], x["START_fw"]) if x["STRAND"] == -1 else (x["START_fw"], x["END_fw"]),
        axis = 1
    ).to_list()
    pfam_domains_genomic_summary = pfam_domains_genomic_summary[['CHROM', 'START', 'END', 'NAME']]
    pfam_domains_genomic_summary["START"] = pfam_domains_genomic_summary["START"].astype(int)
    pfam_domains_genomic_summary["END"] = pfam_domains_genomic_summary["END"].astype(int)

    pfam_domains_genomic_summary.to_csv(output_file, header = False, index = False, sep = '\t')



# ---------------------- CLI Interface ---------------------- #

@click.command()
@click.option('--consensus-panel-rich', type=click.Path(exists=True), help='Input consensus panel file')
@click.option('--domains', type=click.Path(exists=True), help='Input consensus panel file')
@click.option('--output', type=click.Path(), help='Output file')

def main(consensus_panel_rich, domains, output):
    click.echo("Mapping the domains to the genomic sequence...")
    generate_domains2dna_mapping(consensus_panel_rich, domains, output)



if __name__ == '__main__':
    main()


# ---------------------- End of the script ---------------------- #
