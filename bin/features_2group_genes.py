#!/usr/bin/env python


import click
import json
import pandas as pd


def load_table_to_json(table, separator):

    separator_name2character = {"tab": "\t",
                                "comma" : ","
                                }

    # Initialize dictionaries
    link_dict = {}
    gene_dict = {}
    gene_dict_name = {}

    # Read the file
    with open(table, 'r') as file:
        for line in file:
            fields = line.strip().split(separator_name2character.get(separator, "\t"))
            if len(fields) >= 3:
                # Link dictionary: Second field as key and first field as value
                link_dict[fields[1]] = fields[0]

                # Gene dictionary: Second field as key and list of genes as value
                gene_dict[fields[1]] = fields[2:]

                gene_dict_name[fields[0]] = fields[2:]

    return link_dict, gene_dict, gene_dict_name



def load_genes_from_file(file_path):
    '''
    Function to load genes from a file
    '''
    with open(file_path, 'r') as file:
        genes = file.read().splitlines()
    return genes



def filter_genes(data, genes_to_keep_set, minimum_group_size = 1):
    """
    Function to filter genes using set intersection
    """
    filtered_data = {}
    for group, genes in data.items():
        # Use set intersection to filter genes
        filtered_genes = list(set(genes) & genes_to_keep_set)
        if len(filtered_genes) >= minimum_group_size:
            filtered_data[group] = filtered_genes

    return filtered_data


def filtered_groups_to_compressed(data):


    df = pd.DataFrame([(k, v) for k, genes in data.items() for v in [tuple(genes)]], columns=['Pathway', 'Genes'])

    # Create a unique identifier for each set of genes
    df['Genes'] = df['Genes'].apply(tuple)
    df['GeneSetID'] = df['Genes'].apply(lambda x: hash(x))

    # Group pathways by the unique identifier
    gene_set_mapping = df.groupby('GeneSetID').agg({'Pathway' : list, 'Genes': lambda x: set(x).pop()}).reset_index()
    gene_set_mapping["GeneSetID"] = gene_set_mapping["Pathway"].apply(lambda x : "".join([y.capitalize() for y in x[0].split(" ") ]))

    # Create dictionaries from the grouped data
    grouped_pathways = dict(zip(gene_set_mapping['GeneSetID'], gene_set_mapping['Pathway']))
    id_to_genes = dict(zip(gene_set_mapping['GeneSetID'], gene_set_mapping['Genes']))

    groups_to_genes = {g: list(v) for g, v in id_to_genes.items()}

    return groups_to_genes, grouped_pathways





@click.command()
@click.option('--panel-genes-file', required=True, type=click.Path(exists=True), help='Input panel genes file')
@click.option('--add-hotspots', is_flag=True, default=False, help='Flag to add hotspot genes')
@click.option('--hotspot-genes-file', required=False, type=click.Path(exists=True), help='Input hotspot genes file (JSON)')
@click.option('--output-json-groups2names', required=True, type=click.Path(), help='Output JSON file for groups to names')
@click.option('--table-file', required=False, type=click.Path(exists=True), help='Custom gene groups table file')
@click.option('--separator', required=False, type=str, help='Separator for custom gene groups table')
@click.option('--output-json-groups', required=False, type=click.Path(), help='Output JSON file for custom groups')
def main(panel_genes_file, add_hotspots, hotspot_genes_file, output_json_groups2names, table_file, separator, output_json_groups):
    genes_to_keep = set(load_genes_from_file(panel_genes_file))
    genes_to_keep_dict = {x : [x] for x in genes_to_keep}
    genes_to_keep_dict["ALL_GENES"] = [x for x in genes_to_keep if '--' not in x ]

    # handle custom gene groups table
    if table_file and separator and output_json_groups:
        links, gene_groups, gene_groups_name = load_table_to_json(table_file, separator)
        updated_gene_groups = filter_genes(gene_groups_name, genes_to_keep, minimum_group_size= 2)
        grouped_genes_pathways, pathways = filtered_groups_to_compressed(updated_gene_groups)
        genes_to_keep_dict.update(grouped_genes_pathways)
        with open(output_json_groups, 'w') as f:
            json.dump(pathways, f, indent=4)

    # if the hotspots are present add them at the end,
    # this means that if there is any group with the
    # same name as the hotspot it will be overwritten
    if add_hotspots:
        with open(hotspot_genes_file, 'r') as file:
            hotspot_pairs = json.load(file)
            genes_to_keep_dict.update(hotspot_pairs)

    with open(output_json_groups2names, 'w') as f:
        json.dump(genes_to_keep_dict, f, indent=4)

if __name__ == '__main__':
    main()


