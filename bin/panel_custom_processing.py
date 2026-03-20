#!/usr/bin/env python



import click
import polars as pl

from read_utils import custom_na_values_list
muttype_conversion_map = {
                'G/A': 'C/T',
                'G/C': 'C/G',
                'G/T': 'C/A',
                'A/G': 'T/C',
                'A/T': 'T/A',
                'A/C': 'T/G',
            }


def customize_panel_regions(VEP_output_file, custom_regions_file, customized_output_annotation_file,
                            simple = True
                            ):
    """
    Modifies annotations in a VEP output file based on custom genomic regions.

    - For each region in the custom regions file, identifies the corresponding slice
      in the VEP output.
    - Updates gene names and impact values for the region.
    - Saves both the modified annotation file and a record of added regions.

    Args:
        VEP_output_file (str): Path to the full VEP output file (TSV).
        custom_regions_file (str): Custom region definitions (tab-delimited).
        customized_output_annotation_file (str): Output file for updated annotations.
        simple (bool): If True, outputs simplified annotations; else adds more fields.
    """

    # simple = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID'          , 'GENE', 'IMPACT'                                              , 'CONTEXT_MUT', 'CONTEXT']
    # rich   = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID', 'STRAND', 'GENE', 'IMPACT', 'Feature', 'Protein_position', 'Amino_acids', 'CONTEXT_MUT', 'CONTEXT']

    # Read entire file once, partition by chromosome
    all_possible_sites = pl.read_csv(VEP_output_file, separator="\t",
                                        null_values=custom_na_values_list,
                                        schema_overrides={"CHROM": pl.Utf8}
                                    ).with_row_index("__idx")
    print("all possible sites loaded")

    # Get chromosome order before partitioning
    all_chroms = all_possible_sites["CHROM"].unique(maintain_order=True).to_list()

    # Partition by chromosome and release the full dataframe
    chrom_partitions = {
        df["CHROM"][0]: df
        for df in all_possible_sites.partition_by("CHROM", maintain_order=True)
    }
    del all_possible_sites

    custom_regions_df = pl.read_csv(custom_regions_file, separator="\t",
                                        schema_overrides={"CHROM": pl.Utf8})

    # Group custom regions by chromosome
    custom_by_chrom = {}
    for row in custom_regions_df.iter_rows(named=True):
        custom_by_chrom.setdefault(row["CHROM"], []).append(row)

    added_regions_list = []
    write_header = True

    for chrom in all_chroms:
        chr_data = chrom_partitions.pop(chrom)
        print(f"Processing chromosome: {chrom} ({chr_data.height} rows)")

        chr_updates = []

        if chrom in custom_by_chrom:
            for row in custom_by_chrom[chrom]:
                try:
                    # Find start and end indices
                    idx_start = chr_data.filter(pl.col("POS") == row["START"])["__idx"][0]
                    from_start = chr_data.filter(pl.col("__idx") >= idx_start)
                    idx_end = from_start.filter(pl.col("POS") == row["END"])["__idx"][-1]

                    # Extract hotspot data and modify gene names
                    hotspot_data = chr_data.filter(
                        (pl.col("__idx") >= idx_start) & (pl.col("__idx") <= idx_end)
                    ).drop("IMPACT").with_columns(
                        pl.lit(row["NAME"]).alias("GENE")
                    )

                    # Split the string into individual entries
                    entries = row["IMPACTFUL"].split(',')

                    # Create a DataFrame
                    impactful_df = pl.DataFrame({
                        "MUT_ID_pyr": entries,
                        "IMPACT": [row["IMPACT"]] * len(entries)
                    })

                    all_pos_len = hotspot_data.height

                    # convert MUT_ID to C>N and T>N only
                    pyr_expr = pl.col("MUT_ID")
                    for old, new in muttype_conversion_map.items():
                        pyr_expr = pyr_expr.str.replace(old, new, literal=True)
                    hotspot_data = hotspot_data.with_columns(pyr_expr.alias("MUT_ID_pyr"))

                    hotspot_data = hotspot_data.join(impactful_df,
                                                        on = "MUT_ID_pyr",
                                                        how = "full",
                                                        suffix = "_imp"
                                                        )
                    # Coalesce MUT_ID_pyr from both sides of the full join
                    if "MUT_ID_pyr_imp" in hotspot_data.columns:
                        hotspot_data = hotspot_data.with_columns(
                            pl.coalesce(["MUT_ID_pyr", "MUT_ID_pyr_imp"]).alias("MUT_ID_pyr")
                        ).drop("MUT_ID_pyr_imp")

                    all_pos_len_after = hotspot_data.height

                    # TODO add an error raise
                    if all_pos_len != all_pos_len_after:
                        print("Some of the mutations provided are not in the desired region")
                        print(hotspot_data.filter(pl.col("POS").is_null())["MUT_ID_pyr"])
                        hotspot_data = hotspot_data.filter(pl.col("POS").is_not_null())
                        hotspot_data = hotspot_data.with_columns(pl.col("POS").cast(pl.Int64))

                    hotspot_data = hotspot_data.with_columns(
                        pl.col("IMPACT").fill_null(row["NEUTRAL"])
                    ).sort("__idx")

                    ## Collect updates for this chromosome
                    if simple:
                        chr_updates.append(hotspot_data.select(["__idx", "GENE", "IMPACT"]))
                    else:
                        print("Getting Feature to '-'")
                        chr_updates.append(
                            hotspot_data.select(["__idx", "GENE", "IMPACT"]).with_columns(
                                pl.lit("-").alias("Feature")
                            )
                        )

                    added_regions_list.append(hotspot_data)
                    print("Small region added:", row["NAME"])

                except Exception as e:
                    print(f"Error processing row {row}: {e}")

        # Apply updates for this chromosome
        if chr_updates:
            all_updates = pl.concat(chr_updates)
            all_updates = all_updates.unique(subset=["__idx"], keep="last")

            chr_data = chr_data.join(all_updates, on="__idx", how="left", suffix="_upd")

            update_cols = [
                pl.coalesce(["GENE_upd", "GENE"]).alias("GENE"),
                pl.coalesce(["IMPACT_upd", "IMPACT"]).alias("IMPACT"),
            ]
            drop_cols = ["GENE_upd", "IMPACT_upd"]

            if not simple and "Feature_upd" in chr_data.columns:
                update_cols.append(pl.coalesce(["Feature_upd", "Feature"]).alias("Feature"))
                drop_cols.append("Feature_upd")

            chr_data = chr_data.with_columns(update_cols).drop(drop_cols)

        # Deduplicate, sort, and write this chromosome
        chr_data = chr_data.sort("__idx").drop("__idx").unique(
                                                    subset = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID',
                                                                'GENE', 'CONTEXT_MUT', 'CONTEXT', 'IMPACT'],
                                                    keep = 'first',
                                                    maintain_order = True)

        with open(customized_output_annotation_file, 'ab' if not write_header else 'wb') as f:
            chr_data.write_csv(f, separator="\t", include_header=write_header)
        write_header = False

    # Write added regions
    if added_regions_list:
        added_regions_df = pl.concat(added_regions_list, how="diagonal")
        cols_to_drop = [c for c in ["__idx"] if c in added_regions_df.columns]
        if cols_to_drop:
            added_regions_df = added_regions_df.drop(cols_to_drop)
        added_regions_df = added_regions_df.unique(
                                                subset = ['CHROM', 'POS', 'REF', 'ALT', 'MUT_ID',
                                                            'GENE', 'CONTEXT_MUT', 'CONTEXT', 'IMPACT'],
                                                keep = 'first',
                                                maintain_order = True)
        added_regions_df.write_csv('added_regions.tsv' if simple else 'added_regions.rich.tsv',
                                    separator="\t")



@click.command()
@click.option('--vep-output-file', required=True, type=click.Path(exists=True), help='Input VEP output file (TSV)')
@click.option('--custom-regions-file', required=True, type=click.Path(exists=True), help='Input custom regions file (TSV)')
@click.option('--customized-output-annotation-file', required=True, type=click.Path(), help='Output annotation file (TSV)')
@click.option('--simple', is_flag=True, help='Use simple annotation')
def main(vep_output_file, custom_regions_file, customized_output_annotation_file, simple):
    customize_panel_regions(vep_output_file, custom_regions_file, customized_output_annotation_file, simple)

if __name__ == '__main__':
    main()
