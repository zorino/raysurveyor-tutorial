#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# date:    	2016-10-08
# version: 	0.01

import sys
import argparse
import pandas as pd
import numpy as np

from io import StringIO
from Bio import Phylo
from Bio.Phylo import PhyloXML
from Bio.Phylo.TreeConstruction import _DistanceMatrix, _Matrix, DistanceTreeConstructor
from ete3 import Tree, Phyloxml, phyloxml

from scipy.spatial.distance import pdist, squareform

from skbio import DistanceMatrix
from skbio.tree import nj

import plotly.graph_objs as go
import plotly.figure_factory as FF

import colorlover as cl


def read_args():

    parser = argparse.ArgumentParser(description='Create Phylogenetic tree from similarity (Gramian) matrix', add_help=False)

    mandatory_group = parser.add_argument_group(title='Options')
    mandatory_group.add_argument('-i', '--input', metavar='matrix.tsv',
                                 help='Similarity matrix file (Gramian matrix) - Ray-Suveyor output')
    mandatory_group.add_argument('-o', '--output', metavar='tree.nwk', default='tree.nwk', help='Output tree name')
    mandatory_group.add_argument('--dist', default='euclidean',
                                 help='Distance metric [euclidean, cosine, minkowski, .. see https://goo.gl/zg2584')
    mandatory_group.add_argument('--tree', default='nj',
                                 help='Clustering algorithm to build the tree: nj (neighbor joining) or upgma (rooted - not recommended)')
    mandatory_group.add_argument('--norm', help='Normalized the similarity matrix first', action='store_true')
    mandatory_group.add_argument('--norms-matrix', help='Normalized with the norms of this matrix')
    mandatory_group.add_argument('--format', default='newick', help='Output Tree Format [newick or phyloxml]')
    mandatory_group.add_argument('-h', '--help', help='help message', action='store_true')

    args = vars(parser.parse_args())

    if (args['input'] == None) or (args['help']):
        parser.print_help()
        sys.exit()

    return args


def read_matrix(file):
    pd_frame = pd.read_table(file, sep='\t', skipinitialspace=True, index_col=0)
    return pd_frame

def read_categories(file):
    pd_frame = pd.read_table(file, sep='\t', skipinitialspace=True, index_col=None, header=None)
    return pd_frame

def normalize_gram_matrix(K,K2=None):
    normX = None
    if K2 is not None:
        normX = np.sqrt(K2.diagonal())
    else:
        normX = np.sqrt(K.diagonal())
    return ((K/normX).T/normX).T

def condense_matrix(df_dist):
    print("Creating lower triangle of distance matrix..")
    dist_matrix = []
    i = 0
    for index, row in df_dist.iterrows():
        vec = []
        i = i+1
        for j in range(0,i):
            vec.append(float(row[j]))
        dist_matrix.append(vec)

    return dist_matrix

def build_tree(dist_matrix, names_list, clust):

    tree = None
    if clust == 'nj':
        # print(dist_matrix)
        dm = DistanceMatrix(dist_matrix, names_list)
        tree_scikit = nj(dm,result_constructor=str)
        tree = Tree(tree_scikit)
    elif clust == 'upgma':
        dm = _DistanceMatrix(names=names_list, matrix=condense_matrix(dist_matrix))
        constructor = DistanceTreeConstructor()
        tree_biopython = constructor.upgma(dm)
        # remove InnerNode names
        for i in tree_biopython.get_nonterminals():
            i.name = None
        output = StringIO()
        Phylo.write(tree_biopython,output, "newick")
        tree = Tree(output.getvalue())
    else:
        print("Unknown tree clustering method ! Aborting")
        sys.exit()

    return tree

def tree_distances(file):

    t = Tree(file)
    branch_len_out = open(file + ".patristic-dist.tsv", "w")
    avg_distance_leaves = 0

    # Computing patristic distance matrix
    header = ""
    all_leaves = t.get_leaves()
    for i in all_leaves:
        header = header + "\t" + i.name

    nb_of_distances = 0
    max_len = 0
    min_len = 9999999999999999
    branch_len_out.write(header+"\n")
    for leaf1 in all_leaves:
        row = ""
        row += str(leaf1.name)
        for leaf2 in all_leaves:
            distance = np.clip(leaf1.get_distance(leaf2), 0.0, 99999999999999999999999999)
            avg_distance_leaves += distance
            row += "\t%f" % distance
            nb_of_distances += 1
            if distance > max_len:
                max_len = distance
            if distance < min_len and distance > 0:
                min_len = distance

        branch_len_out.write(row+"\n")

    branch_len_out.close()


def plotly_heatmap(gram_matrix, title=None, categories=None, width=600, height=600):

    # Initialize figure by creating upper dendrogram
    figure = None

    figure = FF.create_dendrogram(gram_matrix,
                                  orientation='bottom',
                                  labels=gram_matrix.axes[0])

    for i in range(len(figure['data'])):
        figure['data'][i]['yaxis'] = 'y2'
        figure['data'][i]['showlegend'] = False

    # Create Side Dendrogram
    dendro_side = FF.create_dendrogram(gram_matrix,
                                       orientation='right')

    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
        dendro_side['data'][i]['showlegend'] = False

    # Add Side Dendrogram Data to Figure
    # figure['data'].extend(dendro_side['data'])
    figure.add_traces(dendro_side['data'])

    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))

    data_dist = pdist(gram_matrix)
    heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves,:]
    heat_data = heat_data[:,dendro_leaves]

    heatmap = [
        go.Heatmap(
            x = gram_matrix.axes[0],
            y = gram_matrix.axes[1],
            z = heat_data,
            colorscale = 'RdBu',
            showscale = False,
            # showlegend = False
        )
    ]


    # x_values and labels
    x_values = figure['layout']['xaxis']['tickvals']
    x_labels = figure['layout']['xaxis']['ticktext']
    y_values = dendro_side['layout']['yaxis']['tickvals']

    heatmap[0]['x'] = x_values
    heatmap[0]['y'] = y_values

    # Add Heatmap Data to Figure
    figure.add_traces(heatmap)


    if categories is not None:

        col = categories.columns
        category = categories[col[0]].value_counts()[:11].index.tolist()

        col_p = cl.to_rgb(cl.scales['11']['qual']['Set3'])
        # col_p = cl.to_rgb( cl.scales[str(len(category))]['div']['RdYlBu'] )

        traces = []

        for idx, c in enumerate(category):

            x_categories = categories[categories[col[0]] == c].index.tolist()
            # print(x_categories)
            x_val = [x_values[i] for i, j in enumerate(x_labels) if j in x_categories]

            trace = go.Scatter(
                x = x_val,
                y = [1]*len(x_val),
                mode = "markers",
                marker = {
                    "symbol": "square",
                    "color": col_p[idx],
                    "size": 12
                },
                xaxis = 'x1',
                yaxis = 'y3',
                showlegend = True,
                name = c
            )

            # figure['data'].extend([trace])
            figure.add_trace(trace)


    # Edit Layout
    figure['layout'].update({'width':width, 'height':height, 'title':title,
                             'showlegend': True, 'hovermode': 'closest'})

    # Edit xaxis
    figure['layout']['xaxis'].update({'domain': [.15, 1],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      # "tickmode":False,
                                      "ticks":'',
                                      "tickfont":{
                                          "family":'Arial, sans-serif',
                                          "size": 8,
                                          "color":'black'
                                      }})

    # Edit xaxis2
    figure['layout'].update({'xaxis2': {'domain': [0, .15],
                                        'mirror': False,
                                        'showgrid': False,
                                        'showline': False,
                                        'zeroline': False,
                                        'showticklabels': False,
                                        'ticks':""}})


    # Edit yaxis
    figure['layout']['yaxis'].update({'domain': [0, .85],
                                      'mirror': False,
                                      'showgrid': False,
                                      'showline': False,
                                      'zeroline': False,
                                      'showticklabels': False,
                                      'ticks':""})

    if categories is not None:

        # Edit yaxis2
        figure['layout'].update({'yaxis2':{'domain':[.875, .975],
                                           'mirror': False,
                                           'showgrid': False,
                                           'showline': False,
                                           'zeroline': False,
                                           'showticklabels': False,
                                           'ticks':""}})

        # Edit yaxis3
        figure['layout'].update({'yaxis3':{'domain':[.83, .875],
                                           'mirror': False,
                                           'showgrid': False,
                                           'showline': False,
                                           'zeroline': False,
                                           'showticklabels': False,
                                           'ticks':""}})

    else:

        # Edit yaxis2
        figure['layout'].update({'yaxis2':{'domain':[.825, .975],
                                           'mirror': False,
                                           'showgrid': False,
                                           'showline': False,
                                           'zeroline': False,
                                           'showticklabels': False,
                                           'ticks':""}})


    figure['layout']['yaxis']['ticktext'] = np.asarray(gram_matrix.axes[0][dendro_leaves])
    figure['layout']['yaxis']['tickvals'] = np.asarray(dendro_side['layout']['yaxis']['tickvals'])


    return figure



# Main #
if __name__ == "__main__":

    args = read_args()
    outfile = ""
    outfile = args['output']
    if outfile == "tree.nwk" and args['format'] == 'phyloxml':
        outfile = outfile.replace(".nwk",".xml")

    # read from similarity matrix file
    print("Reading matrix file..")
    df = read_matrix(args['input'])

    print("Computing the distance matrix from gram matrix with [%s] distance metric and Normalized=[%s].."
          % (args['dist'], args['norm']))
    mat_dist = None
    del_minus_fn = lambda *arr: np.clip(*arr,0,99999999999999999999999999)

    if args['norm']:
        norm_dist = None
        if args['norms_matrix']:
            norm_dist = normalize_gram_matrix(df.as_matrix(), read_matrix(args['norms_matrix']).as_matrix())
        else:
            norm_dist = normalize_gram_matrix(df.as_matrix())
        # norm_matrix_out = pd.DataFrame(norm_dist, index=df.index, columns=df.columns)
        # norm_matrix_out.to_csv(path_or_buf=outfile+".normalizedout.tsv", sep="\t")
        mat_dist = squareform(pdist(norm_dist, args['dist']))
        # mat_dist = del_minus_fn(squareform(pdist(norm_dist, args['dist'])))
    else:
        mat_dist = squareform(pdist(df, args['dist']))
        # mat_dist = del_minus_fn(squareform(pdist(df, args['dist'])))

    df_dist = pd.DataFrame(mat_dist, index=df.index, columns=df.columns)

    outfile = args['output']
    df_dist.to_csv(path_or_buf=outfile+".dist.tsv", sep="\t")

    print("Building distance matrix for tree construction..")
    names_list = [str(i).replace('(','_').replace(')','_') for i in df.index.values.tolist()]

    print("Building final [%s] tree.." % args['tree'])
    # tree = build_tree(dist_matrix, names_list, args['tree'])
    tree = build_tree(df_dist, names_list, args['tree'])

    print("Printing tree..")
    if args['format'] == 'phyloxml':
        treedata = tree.write()
        handle = StringIO(treedata)
        tree_out = Phylo.read(handle, "newick")
        phyloxml_out = tree_out.as_phyloxml()
        Phylo.write(phyloxml_out, args['output'], 'phyloxml')
    else:
        tree.write(format=0, outfile=args['output'])
        # Phylo.write(tree, args['output'], 'newick')

    tree_distances(args['output'])
    print("Finish output files : \n  %s\n  %s\n  %s" % (outfile+".dist.tsv", outfile+".patristic-dist.tsv", outfile))
