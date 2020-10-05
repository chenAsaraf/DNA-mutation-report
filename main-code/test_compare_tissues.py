from compare_tissues import compare_tissues
import sys


"""
Test: compare_tissues
"""
def test_compare_tissues():

    healthy = "sample_contigs_k24.contigs.fa"
    tumor = "sample_contigs_k24_tumor.contigs.fa"
    argvlen = len(sys.argv)

    if argvlen > 1:
        test_num = int(sys.argv[1])
    else:
        test_num = 100
        output_prefix = "output"
    if argvlen > 2:
        output_prefix = sys.argv[2]
    else: output_prefix = "output"

    compare_tissues(healthy, tumor, output_prefix, test=True, test_num=test_num)

test_compare_tissues()
