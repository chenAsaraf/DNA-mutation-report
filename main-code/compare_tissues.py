from collections import defaultdict
from Bio import SeqIO
from build_dictionary import tissueDictionary
from matplotlib import pyplot as plt
from mutations_distance import PointMutation
from contigs_analysis import filter_contigs_by_size
import time
import numpy as np

# 1) parse the tumor_cell with non-overlaping window and find there bucket in the dictionary
# 2) search for each contig in the bucket-list in the BST
# 3) find alignment (maby another function)
# 4) send to another function that find the mutation by edit distance

"""
"""
def find_similar_section(tumor_file, output_prefix, k, dictionary, healthyStorage, test=False, test_num=1000):
    # for searching the correct bucket in the dictionary
    # run throw all non-overlaping windows of length k in the sequence (all k-mer)
    filtered_tumor_file, num_tumor_contigs = filter_contigs_by_size(tumor_file, 'filtered_tumor_contigs', test=test, test_num=test_num)
    print("number of filtered tumor contigs:", num_tumor_contigs)
    # initialize the object to save the mutations
    mutations_report = PointMutation(output_prefix)
    records = SeqIO.parse(open(filtered_tumor_file), 'fasta')
    statistics = []
    for tumor_seq in records:
        statistics_of_compares = 0
        contig_len = len(str(tumor_seq.seq))
        for window in range(0, contig_len - k, k): # (last k argument - to jump k chars each step)
            if str(tumor_seq.seq)[window: window + k] in dictionary.keys():
                # Go through all the strings of the healthy tissue in this bucket
                healthy_tissue_records = dictionary[str(tumor_seq.seq)[window: window + k]]
                for record in healthy_tissue_records:
                    healthy_seq = healthyStorage[record.id]
                    # For each alignment - find the overlapping parts and send to the Edit-Distance function
                    for healthy_idx in record.indexes:
                        # Increase the number of compared contigs (For later analysing)
                        statistics_of_compares = statistics_of_compares + 1
                        healthy, tumor = find_overlap(str(healthy_seq), str(tumor_seq.seq), healthy_idx, window)
                        mutations_report.editDistance(tumor, healthy)
        statistics.append(statistics_of_compares)
    avg_compares_per_tumor_seq = sum(statistics)/len(statistics)
    return mutations_report, statistics

def find_overlap(healthy_seq, tumor_seq, healthy_idx, tumor_idx):
    begin_healthy = end_healthy = begin_tumor = end_tumor = 0
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


def plot_distance_histogram(array_of_lengthes, prefix, title="Histogram"):
    fig = plt.figure()
    plt.hist(array_of_lengthes, facecolor='g')
    plt.title(title)
    plt.xlabel('Edit Distance')
    plt.ylabel('Quantity')
    # plt.grid(True)
    fig.savefig(prefix+'.png')



def compare_tissues(healthy_file, tumor_file, output_prefix, test=False, test_num=1000):
    start_time_build_dict = time.time()

    if test:
        dictBuilder = tissueDictionary(healthy_file, test=True, test_num=test_num)
    else:
        dictBuilder = tissueDictionary(healthy_file)

    print()
    print("---------------------------------------")
    print("%.2f seconds to BUILD dictionary " % (time.time() - start_time_build_dict))
    print("---------------------------------------")
    print()
    dictionary, contigsStorage = dictBuilder.get_dictionary_and_storage()
    k = dictBuilder.getK()

    start_time_compare = time.time()
    print()
    print("start to compare tissues...")

    mutations_report, statistics = find_similar_section(tumor_file, output_prefix, k, dictionary, contigsStorage, test=test, test_num=test_num)
    print()
    print("---------------------------------------")
    print("%.2f seconds to COMPARE all tissues " % (time.time() - start_time_compare))
    print("---------------------------------------")
    print()

    sum_comparisons = sum(statistics)
    avg_compares_statistics = sum_comparisons/len(statistics)
    print("number of all comparison:", sum_comparisons)
    print("statistic of contigs comparison per one: ", avg_compares_statistics)

    inserts_amount = 0
    for char, number in mutations_report.inserts.items():
        inserts_amount = inserts_amount + number
    replaces_amount = 0
    for chars, number in mutations_report.replaces.items():
        replaces_amount = replaces_amount + number
    deletes_amount = 0
    for char, number in mutations_report.deletes.items():
        deletes_amount = deletes_amount + number

    print(" ~ Mutations Report ~ ")
    print()
    print("The cancerous tissue from: ", tumor_file)
    print("The noraml tissue from:", healthy_file)
    print()
    print("~ Inserts Amount:", inserts_amount, ", Percentage of all characters: %.2f " % ((float(inserts_amount)/float(mutations_report.sumOfLengths)) * 100), "%")
    print("\t \t", mutations_report.inserts)
    print("~ Replaces Amount:", replaces_amount, ", Percentage of all characters: %.2f " % ((float(replaces_amount)/float(mutations_report.sumOfLengths))*100), "%")
    print("\t \t", mutations_report.replaces)
    print("~ Deletes Amount:", deletes_amount, ", Percentage of all characters: %.2f" % ((float(deletes_amount)/float(mutations_report.sumOfLengths))*100), "%")
    print("\t \t",mutations_report.deletes)
    print("~ Number of comparisons actually entered into the report:", mutations_report.counterOfComparisons)

    # Creating plot of mutation types and save it
    total = plt.figure(figsize=(10, 7))
    keys = inserts.keys()
    values = inserts.values()
    plt.bar(["inserts", "replaces", "deletes"], [inserts_amount, replaces_amount, deletes_amount], align='center', alpha=0.5)
    plt.ylabel('appearance')
    plt.title('Mutation types')
    for i in range(len(list(keys))):
        plt.text(x=i, y=list(values)[i] + 0.1, s=list(values)[i])
    total.savefig(output_prefix + "-total.png")
    
    
#     plt.pie([inserts_amount, replaces_amount, deletes_amount], labels=["inserts", "replaces", "deletes"], autopct='%1.1f%%')
#     inserts = plt.figure(figsize=(10, 7))
#     plt.pie(list(mutations_report.inserts.values()), labels=list(mutations_report.inserts.keys()), autopct='%1.1f%%')
#     replaces = plt.figure(figsize=(10, 7))
#     plt.pie(list(mutations_report.replaces.values()), labels=list(mutations_report.replaces.keys()), autopct='%1.1f%%')
#     deletes = plt.figure(figsize=(10, 7))
#     plt.pie(list(mutations_report.deletes.values()), labels=list(mutations_report.deletes.keys()), autopct='%1.1f%%')
    
#     total.savefig(output_prefix + "-total.png")
#     inserts.savefig(output_prefix + "-inserts.png")
#     replaces.savefig(output_prefix + "-replaces.png")
#     deletes.savefig(output_prefix + "-deletes.png")

    # Create Histogram of all distances for optimizations
    title = "Histogram of all comparisons distances"
    plot_distance_histogram(np.array(mutations_report.listOfDistances), "distances-histogram-" + output_prefix, title)

