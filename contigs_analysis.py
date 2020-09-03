from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def some_function(name, sequence):
    print(name, sequence)

def analayse_histogram(array_of_lengthes):
    bins = np.concatenate(([100,200,300,400,500,1000],np.arange(2000,30000,1000)), axis=None)
    hist, bin_edges = np.histogram(array_of_lengthes, bins = bins)
    fig = plt.figure()
    plt.hist(array_of_lengthes, bins = bins)
    plt.title("Histogram of contgis lengthes")
    #plt.show()
    fig.savefig('histogram_basic24.png')
    return hist, bin_edges

def read_fasta(contigs_file):
    fasta_sequences = SeqIO.parse(open(contigs_file),'fasta')
    counter = 0
    min_length = float('inf')
    max_length = float('-inf')
    sum_length = 0
    list_of_length =[]
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        #new_sequence = some_function(name, sequence)
        # count  the number of contigs
        counter = counter + 1
        contig_len = len(sequence)
        list_of_length.append(contig_len)
        # sum the length of the contigs:
        sum_length = sum_length + contig_len
        # find the minimum and the maximum length of contigs:
        if contig_len < min_length:
            min_length = contig_len
        if contig_len > max_length:
            max_length = contig_len
    # find the avarege length of contigs
    avg = sum_length / counter
    """print("number of contigs:", counter)
    print("min length of contig:", min_length)
    print("max length of contig:", max_length)
    print("avarege length of contigs:", avg)"""
    hist, bin_edges = analayse_histogram(np.array(list_of_length))
    #print(list_of_length)
    """print("histogram:")
    print(hist)
    print("bin_edges:")
    print(bin_edges)"""


read_fasta("../contigs-outputs/basic_k-mer24/basic_try_k-mer24.contigs.fa")
