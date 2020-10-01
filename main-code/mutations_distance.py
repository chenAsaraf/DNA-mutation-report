import edit_distance
from matplotlib import pyplot as plt
import numpy as np

#edit distance between 2 contigs
def editDistance(str1, str2):
    # mistake = len(str1)/10
    sm = edit_distance.SequenceMatcher(a=str1, b=str2)
    # if sm.distance() < mistake:
    f = open("similarContigs.txt", "a")
    f.write(str1 + "\n" + str2 + "\n")
    f.write("The edit distance of this contigs: " + str(sm.distance()) + "\n")
    f.write("The matches of this contigs: " + str(sm.matches()) + "\n")
    temp = [0, 0, 0]
    temp2 = ["The inserts are: ", "The replaces are: ", "The deletes are: "]
    for i in sm.get_opcodes():
        if i[0] == "insert":
            temp[0] += 1
            temp2[0] = temp2[0] + str2[i[3]:i[4]] + ", "
        if i[0] == "replace":
            temp[1] += 1
            temp2[1] = temp2[1] + str2[i[3]:i[4]] + "->" + str1[i[1]:i[2]] + ", "
        if i[0] == "delete":
            temp[2] += 1
            temp2[2] = temp2[2] + str1[i[1]:i[2]] + ", "
    if temp[0] > 0:
        f.write("Number of inserts: " + str(temp[0]) + ". " + str(temp2[0][0:len(temp2[0]) - 2]) + "\n")
    else:
        f.write("Number of inserts: " + str(temp[0]) + ". \n ")
    if temp[1] > 0:
        f.write("Number of replaces: " + str(temp[1]) + ". " + str(temp2[1][0:len(temp2[1])-2]) + "\n")
    else:
        f.write("Number of replaces: " + str(temp[1]) + ". \n ")
    if temp[2] > 0:
        f.write("Number of deletes: " + str(temp[2]) + ". " + str(temp2[2][0:len(temp2[2])-2]) + "\n")
    else:
        f.write("Number of deletes: " + str(temp[2]) + ". \n ")
    f.close()
    editDistance.counters[0] = temp[0] / len(str1)
    editDistance.counters[1] = temp[1] / len(str1)
    editDistance.counters[2] = temp[2] / len(str1)
    editDistance.counters[3] = sm.matches() / len(str1)

editDistance.counters = [0,0,0,0]
str1 = "saturidayi"
str2 = "saturdyj"
editDistance(str1,str2)
print(editDistance.counters)
# Creating plot
fig = plt.figure(figsize=(10, 7))
plt.pie(editDistance.counters, labels=["inserts", "replaces", "deletes", "matches"], autopct='%1.1f%%')

# show plot
plt.show()


