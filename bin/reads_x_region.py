#!/usr/local/bin/python

import re
import sys
import tabix
import multiprocessing
import pandas as pd
import numpy as np

# Define a function to generate exon number
def generate_exon_number(group):
    '''
    TODO
    the way this function is being used will be deprecated
    we should change the implementation
    '''
    group['EXON'] = range(1, len(group) + 1)
    return group


# def parse_mpu(x):
#     '''
#     a typical mpileup output line look like (with option: -s -f $refgen)
#     ..+4ACAC.+4ACAC.+2AC => ['.', '.+ACAC', '.+ACAC', '.+AC']
#     .,,-7tttttgtt => ['.', ',', ',-tttttgt', 't']
#     ,$.,,^>. => [',$', '.', ',', ',', '^.']
#     ,.,.*.**.^*. => [',', '.', ',', '.', '*', '.', '*', '*', '.', '^.']
#     ,....,Tn.t => [',', '.', '.', '.', '.', ',', 'T', 'n', '.', 't']
#     A-1N => ['A-N']
#     '''

#     reads = x["STATUS"].upper()
#     bqlist = x["QNAME"].split(",")

#     readlist = []
#     i = 0        # input pointer in reads
#     j = 0        # output pointer in readlist
#     while i < len(reads):
#         if reads[i] in "ACGTN.,*":
#             readlist.append(reads[i])
#             i += 1
#             j += 1
#         elif reads[i] in '+-':
#             # determine length
#             digit = re.findall('[\+-](\d+)[ACGTN*]+',reads[i:])[0]
#             readlist[j-1] += reads[i] + reads[i+1+len(digit):i+1+len(digit)+int(digit)]
#             i += 1 + len(digit) + int(digit)
#         elif reads[i] == '$':
#             readlist[j-1] += '$'
#             i += 1
#         elif reads[i] == '^':
#             # ^Xa, ^Xa$
#             # readlist.append(reads[i:i+3])
#             readlist.append(reads[i] + reads[i+2])
#             i += 3
#             j += 1
#         else:
#             print('*ERROR* mpileup parser: Unknown char {} in {}[{}]'
#                     .format(reads[i], reads, i), file=sys.stderr)

#     if len(readlist) != len(bqlist):
#         print('*ERROR* mpileup parser: length mismatch between BQ string {} '
#                 'and reads string {} (breakdown={})'
#                 .format(bqlist, reads, ':'.join(readlist)), file=sys.stderr)

#     return readlist, bqlist




def check_if_lost(chromosome, position, old_reads_only, fragments):
    '''
    Check if a given set of reads that are not covering a position any more were lost
    from the previous position or is it that the reads finished
    but the fragment is still there and the mate read is a bit more downstream
    '''
    fragments_of_interest = fragments[fragments["READNAME"].isin(old_reads_only)]
    return sum(np.logical_and(fragments_of_interest["CHROM"] == chromosome, fragments_of_interest["END"] == position))


def compute_in_exon_complete(record_list, frag_data):
    '''
    count unique genomes and unique reads
    '''

    prev_reads = set()
    num_uniq_reads = 0
    uniq_reads = set()
    uniq_genomes = 0

    for record in record_list:
        # Ideally this would be coming from a df or something like this
        # ["CHROM", "POS", "REF", "DEPTH", "STATUS", "QUAL", "QNAME"]
        # print(record[:3])
        # _bps, reads = parse_mpu_single(record[4], record[6])
        # print(reads[:30])
        reads = record[6].split(",")
        chrom = record[0]
        pos = int(record[1])

        set_reads = set(reads)

        new_reads_only = set_reads.difference(uniq_reads)
        old_reads_only = prev_reads.difference(set_reads)

        released_spaces = check_if_lost(chrom, pos, old_reads_only, frag_data)

        len_new = len(new_reads_only)
        num_uniq_reads += len_new

        uniq_reads.update(new_reads_only)

        uniq_genomes += max(0, len_new - (released_spaces) )

        prev_reads = set(reads)
        #print(f'New reads: {len_new}\nReleased genomes: {released_spaces}\nNew genomes: {len_new - released_spaces}')
        #print(f'Total unique genomes: {uniq_genomes}\nTotal unique reads: {num_uniq_reads}\n\n')


    return (uniq_genomes, num_uniq_reads)



def compute_in_exon_complete_df(chrom_pos_reads_values, frag_data):
    '''
    count unique genomes and unique reads
    '''

    prev_reads = set()
    num_uniq_reads = 0
    uniq_reads = set()
    uniq_genomes = 0

    for chrom, pos, set_reads in chrom_pos_reads_values:

        new_reads_only = set_reads.difference(uniq_reads)
        old_reads_only = prev_reads.difference(set_reads)

        released_spaces = check_if_lost(chrom, pos, old_reads_only, frag_data)

        len_new = len(new_reads_only)
        num_uniq_reads += len_new

        uniq_reads.update(new_reads_only)

        uniq_genomes += max(0, len_new - (released_spaces) )

        prev_reads = set_reads
        #print(f'New reads: {len_new}\nReleased genomes: {released_spaces}\nNew genomes: {len_new - released_spaces}')
        #print(f'Total unique genomes: {uniq_genomes}\nTotal unique reads: {num_uniq_reads}\n\n')


    return (uniq_genomes, num_uniq_reads)





# def parse_pileup_rows(list_records):
#     '''
#     Call this function with a list of lists corresponding to lines of a samtools mpileup file
#     and get as output a df with all the rows parsed
#     '''
#     mpileup_data = pd.DataFrame(list_records, columns = ["CHROM", "POS", "REF", "DEPTH", "STATUS", "QUAL", "QNAME"])
#     mpileup_data[["SPLIT_bases", "SPLIT_reads"]] = mpileup_data[["STATUS", "QNAME"]].astype(str).apply(lambda x: pd.Series(parse_mpu(x)), axis = 1)
#     mpileup_data.drop(["STATUS", "QUAL", "QNAME"], axis = 1, inplace = True)
#     return mpileup_data[["CHROM", "POS", "SPLIT_reads"]]




# Function to compute exons info for a subset of regions
def compute_exons_info(sub_regions, tabix_file, frag_data):
    '''
    Wrapper to query the tabix indexed file
    compile all the lines into a single list and
    launch the compute_in_exon_complete_df function
    '''

    chunk_exons_info = {}
    tb = tabix.open(tabix_file)

    for ind, int_row in sub_regions.iterrows():
        query = f"{int_row['CHROM']}:{int_row['START']}-{int_row['END']}"
        try :
            records = tb.querys(query)

            # Collect records into a list
            records_list = []
            for record in records:
                records_list.append( [record[0], int(record[1]), set(record[6].split(",")) ] )

            # Process records to compute exon info
            numgenomes, numreads = compute_in_exon_complete_df(records_list, frag_data)

        except:
            numgenomes = 0
            numreads = 0

        chunk_exons_info[f'{int_row["GENE"]}_{int_row["EXON"]}'] = (numgenomes, numreads)
        print((int_row["GENE"], int_row["EXON"], numgenomes, numreads))

    return chunk_exons_info


# Function to parallelize computation
def parallel_compute(regions, tabix_file, fragments):
    '''
    Wrap up function of the computation of reads for a given panel region
    '''

    num_processes = multiprocessing.cpu_count()  # Adjust as needed
    chunk_size = max(len(regions) // num_processes, 1)
    pool = multiprocessing.Pool(processes=num_processes)
    results = []
    for i in range(num_processes):
        start_idx = i * chunk_size
        end_idx = start_idx + chunk_size if i < num_processes - 1 else len(regions)
        sub_regions = regions.iloc[start_idx:end_idx]
        result = pool.apply_async(compute_exons_info, args=(sub_regions, tabix_file, fragments))
        results.append(result)
    pool.close()
    pool.join()
    exons_info = {}
    for result in results:
        exons_info.update(result.get())
    return exons_info


mpileup_file = sys.argv[1]
fragments_file = sys.argv[2]
regions_file = sys.argv[3]
sample = sys.argv[4]

# ###
# # Read and preprocess the mpileup data
# ###
# mpileup_data = pd.read_csv(mpileup_file, sep = "\t", header = None)
# mpileup_data.columns = ["CHROM", "POS", "REF", "DEPTH", "STATUS", "QUAL", "QNAME"]
# print("mpileup loaded")


# # TODO revise if this can be done more efficiently
# mpileup_data[["SPLIT_bases", "SPLIT_reads"]] = mpileup_data[["STATUS", "QNAME"]].astype(str).apply(lambda x: pd.Series(parse_mpu(x)), axis = 1)
# mpileup_data.drop(["STATUS", "QUAL", "QNAME"], axis = 1, inplace = True)
# print("reads splitted")

regions_df = pd.read_csv(regions_file, sep ='\t', header = 0).iloc[:,:4]
regions_df.columns = ["CHROM", "START", "END", "GENE"]
print("regions loaded")


# Apply the function to each group
regions_df = regions_df.groupby('GENE').apply(generate_exon_number)
print("Exon number generated, numbers based on the coordinate sorted order")

regions_df.to_csv(f"{sample}.custom.bed",
                    sep = "\t",
                    header = True,
                    index = False)

fragments_data = pd.read_csv(fragments_file, sep = "\t", header = None)
fragments_data.columns = ["CHROM", "START", "END", "READNAME"]
print("fragments loaded")


# Call the parallel computation function
exons_info = parallel_compute(regions_df, mpileup_file, fragments_data)




# exons_info = {}
# index = 0
# for ind, row in regions.iterrows():
#     new_end = (mpileup_data.loc[index:,"POS"] > row["END"]).idxmax()
#     numgenomes, numreads = compute_in_exon_complete(mpileup_data.loc[index:new_end+1,["CHROM", "POS", "SPLIT_reads"]])
#     exons_info[f'{row["GENE"]}_{row["EXON"]}'] = (numgenomes, numreads)
#     index = new_end + 1



exons_counts_df = pd.DataFrame.from_dict(exons_info).T
exons_counts_df.columns = ["UNIQ_GENOMES", "UNIQ_READS"]
exons_counts_df = exons_counts_df.reset_index()
exons_counts_df[["GENE", "EXON"]] = exons_counts_df["index"].str.split("_", expand = True)
exons_counts_df = exons_counts_df[['GENE', 'EXON', 'UNIQ_GENOMES', 'UNIQ_READS']]

exons_counts_df_regions = regions_df.merge(exons_counts_df, on = ['GENE', 'EXON'])
exons_counts_df_regions = exons_counts_df_regions[['CHROM', 'START', 'END', 'GENE', 'EXON', 'UNIQ_GENOMES', 'UNIQ_READS']]
exons_counts_df_regions.to_csv(f"{sample}.reads_x_region.tsv.gz",
                                sep = "\t",
                                header = True,
                                index = False)

