from compare_tissues import compare_tissues
import sys


"""
Test: compare_tissues
"""
def test_compare_tissues():

    healthy = "../../contigs-outputs/healthy/basic_k-mer24/basic_try_k-mer24.contigs.fa"
    tumor = "../../contigs-outputs/tumor/basic_k-mer24_T/basic_k-mer24_T.contigs.fa"
    argvlen = len(sys.argv)

    if argvlen > 1:
        output_prefix = sys.argv[1]
    else:
        output_prefix = "output"
        test_num = 100
    if argvlen > 2:
        test_num = sys.argv[2]
    else:
        test_num = None


    compare_tissues(healthy, tumor, output_prefix, test=False, test_num=test_num)

test_compare_tissues()
