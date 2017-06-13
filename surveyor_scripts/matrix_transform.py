#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# date:    	2016-12-15
# version: 	0.01

import sys
import pandas as pd
import numpy as np
import argparse


def read_args():

    parser = argparse.ArgumentParser(description='Make transformation on the Gramian matrix', add_help=False)

    mandatory_group = parser.add_argument_group(title='Options')
    mandatory_group.add_argument('-i', '--input', metavar='matrix.tsv',
                                 help='Similarity matrix file (Gramian matrix) - Ray-Suveyor output')
    mandatory_group.add_argument('-o', '--output', metavar='matrix.transformed.tsv',
                                 default='matrix.transformed.tsv', help='Output matrix name')
    mandatory_group.add_argument('-n', '--normalize', help='Normalize the similarity matrix', action='store_true')
    mandatory_group.add_argument('-m', '--norms-matrix', metavar='norms-matrix.tsv',
                                 help='Use the diagonal of this matrix to normalize the similarity matrix')
    mandatory_group.add_argument('-d', '--drop-indices', metavar='0,1,2', help='Drop indices list separated by coma e.g. 0,1,2')
    mandatory_group.add_argument('-h', '--help', help='help message', action='store_true')

    args = vars(parser.parse_args())

    if (args['input'] == None) or (args['help']):
        parser.print_help()
        sys.exit()

    return args

def read_matrix(file):
    pd_frame = pd.read_table(file,sep='\t', skipinitialspace=True, index_col=0)
    return pd_frame

def drop_indices(matrix, indices):
    new_matrix = matrix.drop(matrix.columns[indices], axis=1)
    new_matrix = new_matrix.drop(matrix.index[indices], axis=0)
    return new_matrix


def normalize_gram_matrix(K,K2=None):
    normX = None
    if K2 is not None:
        normX = np.sqrt(K2.diagonal())
    else:
        normX = np.sqrt(K.diagonal())
    return ((K/normX).T/normX).T


# Main #
if __name__ == "__main__":

    args = read_args()

    outfile = ""
    if args['output']:
        outfile = args['output']
    else:
        outfile = args['input'].replace(".tsv","")

    matrix = read_matrix(args['input'])

    indices = []
    if args['drop_indices']:
        indices = [int(x) for x in args['drop_indices'].split(",")]
        new_matrix = drop_indices(matrix, indices)
    else:
        new_matrix = matrix

    if args['normalize']:
        if args['norms_matrix']:
            if args['drop_indices']:
                norms_matrix = drop_indices(read_matrix(args['norms_matrix']),indices)
            else:
                norms_matrix = read_matrix(args['norms_matrix'])
            norm_matrix = normalize_gram_matrix(new_matrix.as_matrix(),norms_matrix.as_matrix())
        else:
            norm_matrix = normalize_gram_matrix(new_matrix.as_matrix())
    else:
        norm_matrix = new_matrix

    # matrix_new = matrix.drop(matrix.columns[[0]], axis=1)
    # # matrix.drop(matrix.rows[[0,1]], axis=0)
    # matrix_new = matrix_new[1:]

    norm_matrix_out = pd.DataFrame(norm_matrix, index=new_matrix.index, columns=new_matrix.columns)
    norm_matrix_out.to_csv(path_or_buf=outfile, sep="\t")

