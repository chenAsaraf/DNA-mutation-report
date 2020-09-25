from collections import defaultdict
from contigs_analysis import filter_contigs_by_size
from Bio import SeqIO
import sys # for processing toolbar printing


""" Build the dictionary from the healthy cell, for later use """

class tissueDictionary:

    """Inner class- item object to hold the contig information in the dictionary"""
    class dictionaryItem:
        def __init__(self, contig_id, window_index):
            # holds the contig id -
            self.id = int(contig_id)
            # holds list of starting indexes of the k-mer (the key of this entry, the window) in the contig
            self.indexes = [window_index]


    def __init__(self, contigs_file):
        filtered_file = filter_contigs_by_size(contigs_file, 'filtered_'+contigs_file)
        self.k = 10 # window size = 10 for catching 10% of mutations
        # dictionary: holds the k-mer as keys and list of dictionaryItem for each k-mer.
        self.dictionary = defaultdict(list)
        records = SeqIO.parse(open(filtered_file), 'fasta')

        print("start to parsing each contig")
        # setup toolbar
        toolbar_width = 50
        sys.stdout.write("[%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['
        # contigsStorage: holds the contigs sequence in the index of its id

        self.contigsStorage = []
        counter = 0 # this counter is the new id for the contigs. Because maby the contig's file is sample of contigs that shorter then x
        for contig in records:
            self.contigsStorage.append(contig.seq)

            if counter % 100000 == 0:
                sys.stdout.write("-")
                sys.stdout.flush()

            # TODO: delete this line -> print("test: insert the contig number", counter, "to the right cell:", (contig.seq == self.contigsStorage[int(counter)]))
            # parse the windows to dictionary:
            self.parse_window(contig, counter)
            counter = counter + 1

        sys.stdout.write("]\n")  # this ends the progress bar

    def get_dictionary_and_storage(self):
        return self.dictionary, self.contigsStorage

    def getK(self):
        return self.k

    def parse_window(self, contig, counter_id):
    # TODO: make k-mer be a parameter
    # run throw all overlaping windows of length 10 in the sequence (all 10-mer)
        for i in range(len(str(contig.seq))-self.k):
            # key: kmer
            # value: pairs list of (contig id, index of window)
            isExist = False
            # Check if this contig already exist in this entry:
            if str(contig.seq)[i : i+self.k] in self.dictionary.keys():
                entry = self.dictionary[str(contig.seq)[i : i+self.k]]
                # Iterate over all dictionaryItem of this k-mer:
                for value in entry:
                    if value.id == int(counter_id):
                        value.indexes.append(i)
                        isExist = True
                        break
                if not isExist:
                    newItem = self.dictionaryItem(counter_id, i)
                    self.dictionary[str(contig.seq)[i : i+self.k]].append(newItem)
            # Else- if this k-mer is not yet a key:
            else:
                newItem = self.dictionaryItem(counter_id, i)
                self.dictionary[str(contig.seq)[i : i+self.k]].append(newItem)