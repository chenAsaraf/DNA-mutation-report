import edit_distance
import random

#
#
#
class PointMutation:
    def __init__(self):
        self.counters = [0, 0, 0, 0]
        self.counterOfCompares = 0
        self.sumOfLength = 0
        bases = ('A', 'C', 'G', 'T')
        replace = ('AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG')
        self.inserts = dict.fromkeys(bases, 0)
        self.deletes = dict.fromkeys(bases, 0)
        self.replaces = dict.fromkeys(replace, 0)

    def editDistance(self, tumor, healthy):
        mistake = len(tumor)/10  # 10% of mistakes
        print("mistake:" + str(mistake))
        detailsOfDistance = edit_distance.SequenceMatcher(a=tumor, b=healthy)
        print(detailsOfDistance.distance())
        if detailsOfDistance.distance() < mistake:
            self.counterOfCompares += 1
            self.sumOfLength += len(tumor)
            counters = [0, 0, 0]  # counters of: inserts, replaces and deletes
            theChanges = ["The inserts are: ", "The replaces are: ", "The deletes are: "]  # the changes of the Mutations
            for i in detailsOfDistance.get_opcodes():
                if i[0] == "insert":
                    counters[0] += 1
                    theChanges[0] = theChanges[0] + healthy[i[3]:i[4]] + ", "
                    self.inserts[healthy[i[3]:i[4]]] = self.inserts[healthy[i[3]:i[4]]] + 1
                if i[0] == "replace":
                    counters[1] += 1
                    theChanges[1] = theChanges[1] + healthy[i[3]:i[4]] + "->" + tumor[i[1]:i[2]] + ", "
                    temp = healthy[i[3]:i[4]] + tumor[i[1]:i[2]]
                    self.replaces[temp] = self.replaces[temp] + 1
                if i[0] == "delete":
                    counters[2] += 1
                    theChanges[2] = theChanges[2] + tumor[i[1]:i[2]] + ", "
                    temp = str(tumor[i[1]:i[2]])
                    self.deletes[temp] = self.deletes[temp] + 1
            if random.choice(range(1, 101)) > 90:  # take 1/10 from the Mutations to the report
                f = open("similarContigs.txt", "a")
                f.write(tumor + "\n" + healthy + "\n")
                f.write("The edit distance of this contigs: " + str(detailsOfDistance.distance()) + "\n")
                f.write("The matches of this contigs: " + str(detailsOfDistance.matches()) + "\n")
                if counters[0] > 0:
                    f.write("Number of inserts: " + str(counters[0]) + ". " + str(theChanges[0][0:len(theChanges[0]) - 2]) + "\n")
                else:
                    f.write("Number of inserts: " + str(counters[0]) + ". \n")
                if counters[1] > 0:
                    f.write("Number of replaces: " + str(counters[1]) + ". " + str(theChanges[1][0:len(theChanges[1])-2]) + "\n")
                else:
                    f.write("Number of replaces: " + str(counters[1]) + ". \n")
                if counters[2] > 0:
                    f.write("Number of deletes: " + str(counters[2]) + ". " + str(theChanges[2][0:len(theChanges[2])-2]) + "\n")
                else:
                    f.write("Number of deletes: " + str(counters[2]) + ". \n")
                f.close()
            self.counters[0] += counters[0] / len(tumor)
            self.counters[1] += counters[1] / len(tumor)
            self.counters[2] += counters[2] / len(tumor)
            self.counters[3] += detailsOfDistance.matches() / len(tumor)
