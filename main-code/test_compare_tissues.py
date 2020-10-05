from compare_tissues import compare_tissues
import sys


"""
Test: compare_tissues
"""
def test_compare_tissues():
    healthy = "../../contigs-outputs/healthy/basic_k-mer24/basic_try_k-mer24.contigs.fa"
    tumor = "../../contigs-outputs/tumor/basic_k-mer24_T/basic_k-mer24_T.contigs.fa"
    argvlen  = len(sys.argv)
    if argvlen > 1:
        test_num = int(sys.argv[1])
    else:
        test_num = 100
        output_prefix = "output"
    if argvlen > 2:
        output_prefix = sys.argv[2]
    compare_tissues(healthy, tumor, output_prefix, test=True, test_num=test_num)

test_compare_tissues()
