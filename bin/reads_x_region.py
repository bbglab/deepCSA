#!/usr/local/bin/python

import re
import pandas as pd
import numpy as np

import sys

def parse_mpu(x):
    '''
    a typical mpileup output line look like (with option: -s -f $refgen)
    ..+4ACAC.+4ACAC.+2AC => ['.', '.+ACAC', '.+ACAC', '.+AC']
    .,,-7tttttgtt => ['.', ',', ',-tttttgt', 't']
    ,$.,,^>. => [',$', '.', ',', ',', '^.']
    ,.,.*.**.^*. => [',', '.', ',', '.', '*', '.', '*', '*', '.', '^.']
    ,....,Tn.t => [',', '.', '.', '.', '.', ',', 'T', 'n', '.', 't']
    A-1N => ['A-N']
    '''

    reads = x["STATUS"].upper()
    bqlist = x["QNAME"].split(",")

    readlist = []
    i = 0        # input pointer in reads
    j = 0        # output pointer in readlist
    while i < len(reads):
        if reads[i] in "ACGTNacgtn.,*":
            readlist.append(reads[i])
            i += 1
            j += 1
        elif reads[i] in '+-':
            # determine length
            digit = re.findall('[\+-](\d+)[ACGTNacgtn*]+',reads[i:])[0]
            readlist[j-1] += reads[i] + reads[i+1+len(digit):i+1+len(digit)+int(digit)]
            i += 1 + len(digit) + int(digit)
        elif reads[i] == '$':
            readlist[j-1] += '$'
            i += 1
        elif reads[i] == '^':
            # ^Xa, ^Xa$
            # readlist.append(reads[i:i+3])
            readlist.append(reads[i] + reads[i+2])
            i += 3
            j += 1
        else:
            print('*ERROR* mpileup parser: Unknown char {} in {}[{}]'
                    .format(reads[i], reads, i), file=sys.stderr)

    if len(readlist) != len(bqlist):
        print('*ERROR* mpileup parser: length mismatch between BQ string {} '
                'and reads string {} (breakdown={})'
                .format(bqlist, reads, ':'.join(readlist)), file=sys.stderr)

    return readlist, bqlist


def check_if_lost(chromosome, position, old_reads_only, fragments):
    fragments_of_interest = fragments[fragments["READNAME"].isin(old_reads_only)]
    return sum(np.logical_and(fragments_of_interest["CHROM"] == chromosome, fragments_of_interest["END"] == position))


def compute_in_exon_complete(chrom_pos_reads_list):
    prev_reads = set()
    num_uniq_reads = 0
    uniq_reads = set()
    uniq_genomes = 0

    for chrom, pos, reads in chrom_pos_reads_list.values:
        set_reads = set(reads)

        new_reads_only = set_reads.difference(uniq_reads)
        old_reads_only = prev_reads.difference(set_reads)

        released_spaces = check_if_lost(chrom, pos, old_reads_only, fragments_data)

        len_new = len(new_reads_only)
        num_uniq_reads += len_new

        uniq_reads.update(new_reads_only)

        uniq_genomes += max(0, len_new - (released_spaces) )

        prev_reads = set(reads)
        print(f'New reads: {len_new}\nReleased genomes: {released_spaces}\nNew genomes: {len_new - released_spaces}')
        print(f'Total unique genomes: {uniq_genomes}\nTotal unique reads: {num_uniq_reads}\n\n')


    return (uniq_genomes, num_uniq_reads)



mpileup_file = sys.argv[1]
fragments_file = sys.argv[2]
regions_file = sys.argv[3]
sample = sys.argv[4]

###
# Read and preprocess the mpileup data
###
mpileup_data = pd.read_csv(mpileup_file, sep = "\t", header = None)
mpileup_data.columns = ["CHROM", "POS", "REF", "DEPTH", "STATUS", "QUAL", "QNAME"]
print("mpileup loaded")


# TODO revise if this can be done more efficiently
mpileup_data[["SPLIT_bases", "SPLIT_reads"]] = mpileup_data[["STATUS", "QNAME"]].astype(str).apply(lambda x: pd.Series(parse_mpu(x)), axis = 1)
mpileup_data.drop(["STATUS", "QUAL", "QNAME"], axis = 1, inplace = True)
print("reads splitted")

regions = pd.read_csv(regions_file, sep ='\t', header = None).loc[:4,:]
regions.columns = ["CHROM", "START", "END", "GENE"]
print("regions loaded")

# Define a function to generate exon number
def generate_exon_number(group):
    group['EXON'] = range(1, len(group) + 1)
    return group

# Apply the function to each group
regions = regions.groupby('GENE').apply(generate_exon_number)
print("Exon number generated, numbers based on the coordinate sorted order")

fragments_data = pd.read_csv(fragments_file, sep = "\t", header = None)
fragments_data.columns = ["CHROM", "START", "END", "READNAME"]
print("fragments loaded")



exons_info = {}
index = 0
for ind, row in regions.iterrows():
    new_end = (mpileup_data.loc[index:,"POS"] > row["END"]).idxmax()
    numgenomes, numreads = compute_in_exon_complete(mpileup_data.loc[index:new_end+1,["CHROM", "POS", "SPLIT_reads"]])
    exons_info[f'{row["GENE"]}_{row["EXON_ID"]}'] = (numgenomes, numreads)
    index = new_end + 1



exons_counts_df = pd.DataFrame.from_dict(exons_info).T
exons_counts_df.columns = ["UNIQ_GENOMES", "UNIQ_READS"]
exons_counts_df = exons_counts_df.reset_index()
exons_counts_df[["GENE", "EXON"]] = exons_counts_df["index"].str.split("_", expand = True)
exons_counts_df = exons_counts_df[['GENE', 'EXON', 'UNIQ_GENOMES', 'UNIQ_READS']]
exons_counts_df.to_csv(f"{sample}.reads_x_region.tsv.gz",
                            sep = "\t",
                            header = True,
                            index = False)

