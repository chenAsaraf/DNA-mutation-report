from collections import defaultdict
from Bio import SeqIO
from build_dictionary import tissueDictionary

# 1) parse the tumor_cell with non-overlaping window and find there bucket in the dictionary
# 2) search for each contig in the bucket-list in the BST
# 3) find alignment (maby another function)
# 4) send to another function that find the mutation by edit distance

"""
"""
def parse_tumor_with_window(tumor_file, k, dictionary, healthyStorage):
    # for searching the correct bucket in the dictionary
    # run throw all NON-overlaping windows of length k in the sequence (all k-mer)
    parsed_tumor_file = SeqIO.parse(open(tumor_file), 'fasta')
    for tumor_seq in parsed_tumor_file:
        contig_len = len(str(tumor_seq.seq))
        for window in range(0, contig_len - k, k): # (last k argument - to jump k chars each step)
            if str(tumor_seq.seq)[window: window + k] in dictionary.keys():
                # Go through all the strings of the healthy tissue in this bucket
                healthy_tissue_records = dictionary[str(tumor_seq.seq)[window: window + k]]
                for record in healthy_tissue_records:
                    healthy_seq = healthyStorage[record.id]
                    # For each alignment - find the overlapping parts and send to the Edit-Distance function
                    for healthy_idx in record.indexes:
                        healthy, tumor = find_overlap(healthy_seq, tumor_seq, healthy_idx, window)


def find_overlap(healthy_seq, tumor_seq, healthy_idx, tumor_idx):
    begin_healthy, end_healthy, begin_tumor, end_tumor = 0
    # 1) Find the beginning of the overlap:
    # If one or both of the indexes is 0 - the start of the sequence, leave the indexes the same
    if tumor_idx == 0 or healthy_idx == 0:
        begin_healthy = healthy_idx
        begin_tumor = tumor_idx
    # Else- if one is larger then the other - take the sequence with the smaller index
    # from the beginning, and take the sequence with the larger index from:
    # larger_index - smaller_index
    else:
        if healthy_idx < tumor_idx:
            begin_healthy = 0
            begin_tumor = tumor_idx - healthy_idx
        else:
            begin_tumor = 0
            begin_healthy = healthy_idx - tumor_idx
    # 2) Find the ending of the overlap:
    remaining_healthy = len(healthy_seq) - healthy_idx
    remaining_tumor = len(tumor_seq) - tumor_idx
    if remaining_healthy < remaining_tumor:
        end_healthy = healthy_idx + remaining_healthy
        end_tumor = tumor_idx + remaining_healthy
    else:
        end_healthy = healthy_idx + remaining_tumor
        end_tumor = tumor_idx + remaining_tumor
    return healthy_seq[begin_healthy : end_healthy], tumor_seq[begin_tumor : end_tumor]


def compare_tissues(healthy_file, tumor_file):
    dictBuilder = tissueDictionary(healthy_file)
    dictionary, contigsStorage = dictBuilder.get_dictionary_and_storage()
    k = dictBuilder.getK()
    parse_tumor_with_window(tumor_file, k, dictionary, contigsStorage)