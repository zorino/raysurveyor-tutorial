#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# date:    	2016-10-18
# version: 	0.01

# import to compare trees
import sys
import argparse
from ete3 import Tree

# import to compare clusters
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.stats import pearsonr
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd

from skbio import DistanceMatrix
from skbio.stats.distance import mantel


def read_args():

    parser = argparse.ArgumentParser(description='Compare phylogenetic trees and clusters drawn from distance matrices',
                                     add_help=False)

    mandatory_group = parser.add_argument_group(title='Options')
    mandatory_group.add_argument('-t', '--type', default='trees',
                                 help='trees or clusters (default=trees)')
    mandatory_group.add_argument('-o', '--output', metavar='tree-compare.tsv',
                                 default='tree-compare.tsv', help='Output comparison file')
    mandatory_group.add_argument('-r', '--ref', metavar='ref tree/matrix',
                                 help='Reference tree or dist matrix')
    mandatory_group.add_argument('-a', '--all', metavar='trees/matrices...', action='append',
                                 help='All other trees or dist matrices to compare with')
    mandatory_group.add_argument('-h', '--help', help='help message', action='store_true')

    args = vars(parser.parse_args())

    if (args['type'] == None or
        args['ref'] == None or
        args['all'] == None or
        args['help']):
        parser.print_help()
        sys.exit()

    return args


def compare_trees(args):

    tree_ref = Tree(args['ref'])
    trees = []
    for i in args['all']:
        tree = {}
        tree['name'] = i
        tree['comparison'] = tree_ref.compare(Tree(i),unrooted=True)
        trees.append(tree)

    with open(args['output'],"w") as f_out:
        f_out.write("Tree\tEffective_Tree_Size\tRF-normalized\tRF\tRF-max\t% edges in Src\t% edges in Ref\n")
        for tree in trees:
            f_out.write(tree['name']+"\t")
            f_out.write(str(tree['comparison']['effective_tree_size'])+"\t")
            f_out.write(str(tree['comparison']['norm_rf'])+"\t")
            f_out.write(str(tree['comparison']['rf'])+"\t")
            f_out.write(str(tree['comparison']['max_rf'])+"\t")
            f_out.write(str(tree['comparison']['ref_edges_in_source'])+"\t")
            f_out.write(str(tree['comparison']['source_edges_in_ref'])+"\t")
            f_out.write("\n")

    # print(trees[0]['comparison'])


def check_symmetry(m):
    # will also correct different number
    for i in range(0,len(m)-1):
        for j in range(0,len(m)-1):
            try:
                val = float(m[i][j])
                val2 = float(m[j][i])
            except ValueError:
                print("Not an number!")
            if(m[i][j] != m[j][i]):
                m[j][i] = m[i][j]


def mantel_test(ref, query):

    return mantel(ref, query)


def cophentic_correlation_test(ref, query):

    ref = ref.as_matrix()
    query = query.as_matrix()
    linkage_ref = linkage(pdist(ref), 'average')
    c_ref, coph_dists_ref = cophenet(linkage_ref, pdist(ref))
    linkage_fst = linkage(pdist(query), 'average')
    c_fst, coph_dists_query = cophenet(linkage_fst, pdist(query))
    cophenetic_pearson, p_value_cophenetic = pearsonr(coph_dists_ref, coph_dists_query)

    return cophenetic_pearson, p_value_cophenetic

def compare_clusters(args):

    ref_df = pd.read_table(args['ref'], sep='\t', skipinitialspace=True, index_col=0).as_matrix()
    check_symmetry(ref_df)
    linkage_ref = linkage(ref_df, 'average')
    c_ref, coph_dists_ref = cophenet(linkage_ref, pdist(ref_df))

    outfile = open(args['output'],"w")
    outfile.write("Tree_cluster\tMantel_Correlation_Coefficient\tManter_P-value\tCophenetic_Pearson\tCophenetic_P-value\n")

    for i in args['all']:
        fst_df = pd.read_table(i, sep='\t', skipinitialspace=True, index_col=0).as_matrix()
        check_symmetry(fst_df)
        mantel_coeff = 0.0
        p_value_mantel = 0.0
        cophenetic_pearson = 0.0
        p_value_cophenetic = 0.0
        n = 0
        try:
            # mantel_coeff, p_value_mantel, n = mantel(ref_df, fst_df)
            mantel_coeff, p_value_mantel, n = mantel_test(ref_df, fst_df)
            linkage_fst = linkage(fst_df, 'average')
            c_fst, coph_dists_fst = cophenet(linkage_fst, pdist(fst_df))
            cophenetic_pearson, p_value_cophenetic = pearsonr(coph_dists_ref, coph_dists_fst)
        except Exception as e:
            print("Error : %s" % str(e))
            mantel_coeff = "Failed"
            p_value_manel = "Failed"
            cophenetic_pearson = "Failed"
            p_value_cophenetic = "Failed"

        outfile.write(i+"\t"+str(mantel_coeff)+"\t"+str(p_value_mantel)+"\t"+str(cophenetic_pearson)+"\t"+str(p_value_cophenetic)+"\n")

    outfile.close()


# Main #
if __name__ == "__main__":

    args = read_args()

    if args['type'] == "trees":
        compare_trees(args)
    elif args['type'] == "clusters":
        compare_clusters(args)
    else:
        print("Bad Type ! .. trees or clusters")
        sys.exit()

    # print("Compare tree")
    # print(tree_ref.compare(trees[0],unrooted=True))
    # # print(trees[0].name)

