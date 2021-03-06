from compare_tissues import find_overlap

""" Test : find_overlap:
 case 1:  If one or both of the indexes is 0
    window = "GGGCTA", k = 6
    string A = "ACTCTGGGCTAACTTTTT", index = 5
    string B =      "GGGCTAACTTTTTGGCTACC", index = 0
    results should be: "GGGCTAACTTTTT" x 2
 """
def test_1_find_overlap():
    resultA, resultB = find_overlap("ACTCTGGGCTAACTTTTT", "GGGCTAACTTTTTGGCTACC", 5, 0)
    print("Expected results: GGGCTAACTTTTT, The actual results: 1-", resultA, "2-", resultB)
    print("The correct results?", resultA == "GGGCTAACTTTTT")

"""
 case 2:  None of the first indexes is 0
    window = "GGGCTAACTT", k = 10
    string A = "ACTCTGGGCTAACTTTTT", index = 5
    string B =   "TCTGGGCTAACTTTTTGGCTACC", index = 3
    results should be: "TCTGGGCTAACTTTTT" x 2
 """
def test_2_find_overlap():
    resultA, resultB = find_overlap("ACTCTGGGCTAACTTTTT", "TCTGGGCTAACTTTTTGGCTACC", 5, 3)
    print("Expected results: TCTGGGCTAACTTTTT, The actual results: 1-", resultA, "2-", resultB)
    print("The correct results?", resultA == "TCTGGGCTAACTTTTT")


"""
 case 3: One of the strings is contained in the other
    window = "TGGG", k = 4
    string A = "ACTCTGGGCTAACTTTTT", index = 4
    string B =  "CTCTGGGCTAA", index = 3
    results should be: "CTCTGGGCTAA" x 2
 """
def test_3_find_overlap():
    resultA, resultB = find_overlap("ACTCTGGGCTAACTTTTT", "CTCTGGGCTAA", 4, 3)
    print("Expected results: CTCTGGGCTAA, The actual results: 1-", resultA, "2-", resultB)
    print("The correct results?", resultA == "CTCTGGGCTAA")



test_1_find_overlap()
print("------------")
test_2_find_overlap()
print("------------")
test_3_find_overlap()
