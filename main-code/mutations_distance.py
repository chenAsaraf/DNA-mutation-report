import edit_distance
import random
import math
from Bio import pairwise2


INSERTS = 0
REPLACES = 1
DELETES = 2


class PointMutation:
    """
    This class is based on the 'edit_distance' library, that implements
    Levenshtein distance.

    Informally, the Levenshtein distance between two words is the minimum
    number of single-character edits (insertions, deletions or
    substitutions) required to change one word into the other (by wikipedia)

   We used the following documentation:
   https://docs.python.org/2/library/difflib.html

    Attributes
    ----------
    output_prefix : string
        represents the output file prefix
    counterOfComparisons : int
        variable to count the number of all compared contigs
    sumOfLengths : int
        variable to sum the lengths of all the compared contigs
    listOfDistances : list(int)
        list of all the edit-distance comparisons results
    inserts : dictionary
        each key in the dictionary is a char from the set {'A', 'C', 'G', 'T'}
        each value is the amount of characters of that type that have inserted
    deletes : dictionary
        each key in the dictionary is a char from the set {'A', 'C', 'G', 'T'}
        each value is the amount of characters of that type that have deleted
    replaces : dictionary
        each key in the dictionary is a char from the set
        {'AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG'}
        the first char in each key represents the replaced char
        the second char in each key represents the char inserted in place
        i.e 'AC' represents that 'A' was in the healthy tissue and was replaced
        with 'C' in the tumor tissue each value is the amount of characters of
        that type that have replaces
    in_percentages : float
        the percentages of chars that inserted out of all the chars we tested
    del_percentages : float
        the percentages of chars that deleted out of all the chars we tested
    rep_percentages : float
        the percentages of chars that replaced out of all the chars we tested

    Methods
    -------
    editDistance(healthy, tumor)

    set_in_percentages(perc)

    set_del_percentages(perc)

    set_rep_percentages(perc)

    """

    def __init__(self, output_prefix):
        """

        :param output_prefix: string
                    represents the output file prefix
        """
        self.output_prefix = output_prefix
        self.counterOfComparisons = 0
        self.sumOfLengths = 0
        self.listOfDistances = [] #
        # The keys for the Dictionary of inserts and deletes:
        bases = ('A', 'C', 'G', 'T')
        # The keys for the Dictionary of replaces
        replace = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')
        # Set all the values in the Dictionary to 0
        self.inserts = dict.fromkeys(bases, 0)
        self.deletes = dict.fromkeys(bases, 0)
        self.replaces = dict.fromkeys(replace, 0)
        self.in_percentages = 0
        self.del_percentages = 0
        self.rep_percentages = 0

    def tissues_edit_distance(self, healthy, tumor):
        """Calculate Levenshtein distance to the healthy and tumor strings

        :param healthy: string
                represents the genome from the healthy tissue
        :param tumor: string
                represents the genome from the tumor tissue
        """
        # We assume less than 15 percent of the distance for being the same tissue:
        errors_percent = math.ceil(len(tumor) * 0.15)
        # The first string (i.e 'a') is the one that is compared to the second.
        # The second string (i.e 'b') does not change.
        sequence_match = edit_distance.SequenceMatcher(a=healthy, b=tumor)
        distance = sequence_match.distance()  # the edit distance between the strings
        self.listOfDistances.append(distance)
        if distance < errors_percent:
            self.counterOfComparisons += 1
            self.sumOfLengths += len(tumor)
            # Counters of: inserts, replaces and deletes
            counters_for_mutations = [0, 0, 0]
            # Strings for later printed samples:
            the_inserts = "Actual inserts: "
            the_replaces = "Actual replaces: "
            the_deletes = "Actual deletes: "
            # get_opcodes() returns list of 5-tuples describing how to turn a into b
            for opcode in sequence_match.get_opcodes():
                if opcode[0] == "insert":
                    counters_for_mutations[INSERTS] += 1
                    inserted_char = tumor[opcode[3]:opcode[4]]
                    # Add the char that inserted to the tumor string
                    the_inserts = the_inserts + inserted_char + ", "
                    # Increase the char-entry in the dictionary
                    self.inserts[inserted_char] += 1
                if opcode[0] == "replace":
                    counters_for_mutations[REPLACES] += 1
                    replaced_from = healthy[opcode[1]:opcode[2]]
                    replaced_to = tumor[opcode[3]:opcode[4]]
                    the_replaces = the_replaces + replaced_from + "->" + replaced_to + ", "
                    # Increase the chars-entry in the dictionary
                    self.replaces[replaced_from + replaced_to] += 1
                if opcode[0] == "delete":
                    counters_for_mutations[DELETES] += 1
                    deleted_char = healthy[opcode[1]:opcode[2]]
                    the_deletes = the_deletes + deleted_char + ", "
                    # Increase the chars-entry in the dictionary
                    self.deletes[deleted_char] += 1
            # One in 100 strings enters a report that prints the comparable sections themselves:
            if random.choice(range(1, 101)) > 90 and distance > 0:
                f = open(self.output_prefix + ".txt", "a")  # The sampling report
                f.write("\n============================================================================\n")
                # Print the comparable sections with alignment:
                alignment = pairwise2.align.globalxx(tumor, healthy)
                f.write("The Tumor Tissue: \n" + alignment[0].seqA + "\nThe Healthy Tissue:\n" + alignment[0].seqB + "\n")  # Add the contigs to the file
                # Print the mutations of each type:
                if counters_for_mutations[INSERTS] > 0:
                    f.write("Inserts Amount: " + str(counters_for_mutations[INSERTS]) + ". " + the_inserts[0:len(the_inserts)-2] + "\n")
                else:
                    f.write("Inserts Amount: " + str(counters_for_mutations[INSERTS]) + ". \n")
                if counters_for_mutations[REPLACES] > 0:
                    f.write("Replaces Amount: " + str(counters_for_mutations[REPLACES]) + ". " + the_replaces[0:len(the_replaces)-2] + "\n")
                else:
                    f.write("Replaces Amount: " + str(counters_for_mutations[REPLACES]) + ". \n")
                if counters_for_mutations[DELETES] > 0:
                    f.write("Deletes Amount: " + str(counters_for_mutations[DELETES]) + ". " + the_deletes[0:len(the_deletes)-2] + "\n")
                else:
                    f.write("Deletes Amount: " + str(counters_for_mutations[DELETES]) + ". \n")
                f.close()

    def set_in_percentages(self, perc):
        self.in_percentages = perc

    def set_del_percentages(self, perc):
        self.del_percentages = perc

    def set_rep_percentages(self, perc):
        self.rep_percentages = perc
