from collections import defaultdict
from Bio import SeqIO
import sys # for processing toolbar printing
from contigs_analysis import filter_contigs_by_size


class TissueDictionaryBuilder:
    """
    This class builds a dictionary from the healthy tissue


    Construction process:

    1. filter_contigs _by _size -
       We select the contigs from the healthy tissue in a certain length range as
       an initial implementation and based on a histogram of lengths that we
       extracted in advance from the file analysis. This should be improved in
       future versions to be capable analyse a wider range of lengths.

    2. For each contig in the filtered file:

        2.1. Append the contig to the contigsStorage -
             we use Python's List to store all the sequences of the contigs.
             The index of the contig in the list is its ID, and thus later in the
             program it will allow quick access to each contig by the "[]" operator.

        2.2. __parse_window -
             Window parsing of the contig sequence: A moving window of 10 chars
             move along the sequence, inserts the index of the sequence into the
             dictionary. The content of the window become the key and the value is id.


    Attributes
    ----------
    k : int
        represent the window size (the kmer size). equal to 10
    dictionary : defaultdict(list)
        holds the k-mer as keys and list of dictionaryItem as values
    contigsStorage : list
        holds the contigs sequence in the index of its id


    Methods
    -------
    get_dictionary_and_storage()
        :return self.dictionary, self.contigsStorage

    get_k()
    :return self.k

    """

    class DictionaryItem:
        """
        Inner class- item object to hold the contig information in the dictionary

        Attributes
        ----------
        id : int
            holds the contig id
        indexes : list(int)
            holds list of starting indexes of the k-mer
            (the key of this entry, the window) in the contig
        """
        def __init__(self, contig_id, window_index):
            """
            :param contig_id: int
            :param window_index: int
            """
            self.id = int(contig_id)
            self.indexes = [window_index]


    def __init__(self, contigs_file, test=False, test_num=1000):
        """

        :param contigs_file: 'FASTA' file
                    healthy tissue contigs file from which to build
                    the dictionary
        :param test: boolean (optional, default = False)
                    a variable designed to assist in the software
                    development process. If test=True then the
                    software will only run up to test_num contigs.
        :param test_num: int (optional, default = 1000)
                    this parameter used only in case test=True
        """
        filtered_file, num_contigs = filter_contigs_by_size(contigs_file, 'filtered_contigs', test=test, test_num=test_num)
        print("number of filtered contigs:", num_contigs)
        self.k = 10
        self.dictionary = defaultdict(list)
        records = SeqIO.parse(open(filtered_file), 'fasta')

        print("start to parsing each contig")
        # setup toolbar of process progress
        if test_num > 50:
            toolbar_width = 50
            eval_progress = int(num_contigs / 50)
        if test_num < 51:
            toolbar_width = 10
            eval_progress = 2
        sys.stdout.write("[%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['

        self.contigsStorage = []
        # The following counter is the new ID for each contig.
        # This is because the contig ID in the filtered file may not be
        # in a consistent order. We want the ID order to be consistent
        # for later on access the sequence in the storage list by i'ts ID index.
        counter = 0
        for contig in records:
            self.contigsStorage.append(contig.seq)

            if counter % eval_progress == 0:
                sys.stdout.write("-")
                sys.stdout.flush()

            # parse the windows to dictionary:
            self.__parse_window(contig, counter)
            counter = counter + 1

        sys.stdout.write("]\n")  # this ends the progress bar

    def __parse_window(self, contig, counter_id):
        """ Move the window along the contig sequence and enter the
         contig information at the appropriate entries in the dictionary


        :param contig: SeqRecord object
            biopython object - hold the contig sequence and identifiers
        :param counter_id: int
        """
        # Run throw all overlapping windows of length 10 in the sequence (all 10-mer)
        for i in range(len(str(contig.seq))-self.k):
            # Check if this contig already exist in this entry:
            is_exist = False
            if str(contig.seq)[i : i+self.k] in self.dictionary.keys():
                entry = self.dictionary[str(contig.seq)[i : i+self.k]]
                # Iterate over all DictionaryItem of this k-mer:
                for value in entry:
                    if value.id == int(counter_id):
                        value.indexes.append(i)
                        is_exist = True
                        break
                if not is_exist:
                    new_item = self.DictionaryItem(counter_id, i)
                    self.dictionary[str(contig.seq)[i : i+self.k]].append(new_item)
            # Else- if this k-mer is not yet a key:
            else:
                new_item = self.DictionaryItem(counter_id, i)
                self.dictionary[str(contig.seq)[i : i+self.k]].append(new_item)

    def get_dictionary_and_storage(self):
        return self.dictionary, self.contigsStorage

    def get_k(self):
        return self.k

