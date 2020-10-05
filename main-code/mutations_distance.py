import edit_distance
import random
import math
from Bio import pairwise2


"""
This modulo is based on the 'edit_distance' library, that implements the Levenshtein distance.
Informally, the Levenshtein distance between two words is the minimum number of single-character edits
(insertions, deletions or substitutions) required to change one word into the other (wikipedia)
To be noted: The 'edit_distance' library implementation could likely be optimized to be faster
within Python. And could probably be much faster if implemented in C.
"""

class PointMutation:
    def __init__(self, output_prefix):
        self.output_prefix = output_prefix
        self.counters = [0, 0, 0, 0]  # counters for each mutation types: inserts, replaces, deletes and matches
        self.counterOfCompares = 0  # number of compared contigs
        self.sumOfLength = 0  # sum length of all the compared contigs (for avg)
        bases = ('A', 'C', 'G', 'T')  # the keys for the Dictionary of inserts and deletes
        replace = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')  # the keys for the Dictionary of replaces
        # set all the values in the Dictionary to 0
        self.inserts = dict.fromkeys(bases, 0)
        self.deletes = dict.fromkeys(bases, 0)
        self.replaces = dict.fromkeys(replace, 0)

    def editDistance(self, tumor, healthy):
        errors_precent = math.ceil(len(tumor)/10)  # 10% of mistakes
        sequence_match = edit_distance.SequenceMatcher(a=tumor, b=healthy)
        distance = sequence_match.distance()  # distance between the contigs
        matches = sequence_match.matches()  # the matches (the same chars) between the contigs

        if distance < errors_precent:
            self.counterOfCompares += 1
            self.sumOfLength += len(tumor)
            counters = [0, 0, 0]  # counters of: inserts, replaces and deletes
            theChanges = ["The inserts are: ", "The replaces are: ", "The deletes are: "]  # the changes of the Mutations
            for opcode in sequence_match.get_opcodes():  # Return list of 5-tuples describing how to turn a into b
                if opcode[0] == "insert":
                    counters[0] += 1
                    theChanges[0] = theChanges[0] + healthy[opcode[3]:opcode[4]] + ", "  # Add the char that we insert to the tumor contig
                    self.inserts[healthy[opcode[3]:opcode[4]]] += 1  # Add 1 to the values of this key
                if opcode[0] == "replace":
                    counters[1] += 1
                    theChanges[1] = theChanges[1] + healthy[opcode[3]:opcode[4]] + "->" + tumor[opcode[1]:opcode[2]] + ", "
                    self.replaces[healthy[opcode[3]:opcode[4]] + tumor[opcode[1]:opcode[2]]] += 1  # Add 1 to the values of this key
                if opcode[0] == "delete":
                    counters[2] += 1
                    theChanges[2] = theChanges[2] + tumor[opcode[1]:opcode[2]] + ", "
                    self.deletes[tumor[opcode[1]:opcode[2]]] += 1  # Add 1 to the values of this key
            if random.choice(range(1, 101)) > 90 and distance > 0:  # Take 1/10 from the compared contigs to the sampling report in term that there are at least 1 mutation between them
                f = open(self.output_prefix + ".txt", "a")  # The sampling report
                f.write("\n============================================================================\n")
                alignment = pairwise2.align.globalxx(tumor, healthy)
                f.write(pairwise2.format_alignment(*alignment[0]))
                # f.write("Edit Distance: " + str(distance) + "\n")
                # f.write("The matches of this contigs: " + str(sequence_match.matches()) + "\n")
                # For each mutation type, check if there are mutations like this type and add to the sampling report with the mutations themselves
                if counters[0] > 0:
                    f.write("Inserts: " + str(counters[0]) + ". " + str(theChanges[0][0:len(theChanges[0]) - 2]) + "\n")
                else:
                    f.write("Inserts: " + str(counters[0]) + ". \n")
                if counters[1] > 0:
                    f.write("Replaces: " + str(counters[1]) + ". " + str(theChanges[1][0:len(theChanges[1])-2]) + "\n")
                else:
                    f.write("Replaces: " + str(counters[1]) + ". \n")
                if counters[2] > 0:
                    f.write("Deletes: " + str(counters[2]) + ". " + str(theChanges[2][0:len(theChanges[2])-2]) + "\n")
                else:
                    f.write("Deletes: " + str(counters[2]) + ". \n")
                f.close()
            self.counters[0] += counters[0] / len(tumor)  # Add the part of inserts in the contig
            self.counters[1] += counters[1] / len(tumor)  # Add the part of replaces in the contig
            self.counters[2] += counters[2] / len(tumor)  # Add the part of deletes from the contig
            self.counters[3] += matches / len(tumor)  # Add the part of matches between the contigs
