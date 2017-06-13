#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# date:    	2016-10-05
# version: 	0.01


import sys
import os
from subprocess import call
import pyvolve
import pyani
from ete3 import Tree
from Bio import SeqIO


def split_seq(seq_f,outdir):
    os.mkdir(outdir) if not os.path.exists(outdir) else None
    handle = open(seq_f, "r")
    for record in SeqIO.parse(handle, "fasta"):
        output_handle = open(outdir+"/"+record.id.replace("/","_")+".fasta", "w")
        SeqIO.write(record,output_handle,"fasta")
        output_handle.close()
    handle.close()

def pyani_seq(seq_f):
    outdir = seq_f.replace(".fasta","")
    print(outdir)
    split_seq(seq_f,outdir)
    call("average_nucleotide_identity.py --workers 4 -i %s -o %s -m ANIb" % (outdir, outdir+".ANI"), shell=True)

def tree_distances_info(file,scale,seq_len):

    t = Tree(file)
    # branch_len_matrix_f = file + ".branches-len.tsv"
    branch_len_out = open(file + ".%d.patristic-dist.tsv" % seq_len, "w")
    tree_info = open(file + ".%d.tree-info.txt" % seq_len, "w")
    avg_distance_leaves = 0

    # Computing patristic distance matrix
    header = ""
    all_leaves = t.get_leaves()
    for i in all_leaves:
        header = header + "\t" + i.name

    nb_of_distances = 0
    max_len = 0
    min_len = 99999999999999
    branch_len_out.write(header+"\n")
    for leaf1 in all_leaves:
        row = ""
        row += str(leaf1.name)
        for leaf2 in all_leaves:
            avg_distance_leaves += leaf1.get_distance(leaf2)
            distance = leaf1.get_distance(leaf2)
            row += "\t%f" % distance
            nb_of_distances += 1
            if distance > max_len:
                max_len = distance
            if distance < min_len and distance > 0:
                min_len = distance

        branch_len_out.write(row+"\n")

    tree_info.write("Scale_factor(1=original-tree)\t%f\n" % scale)
    tree_info.write("Seq_Length\t%d\n" % seq_len)
    tree_info.write("Number_of_leaves_(taxa)\t%d\n" % len(all_leaves))
    tree_info.write("Minimal_patristic_distance\t%f\n" % min_len)
    tree_info.write("Maximal_patristic_distance\t%f\n" % max_len)
    tree_info.write("Average_patristic_distance\t%f\n" % (avg_distance_leaves/(nb_of_distances*scale)))

    print("Scale_factor(1=original-tree)\t%f" % scale)
    print("Seq_Length\t%d" % seq_len)
    print("Number_of_leaves_(taxa)\t%d" % len(all_leaves))
    print("Minimal_patristic_distance\t%f" % min_len)
    print("Maximal_patristic_distance\t%f" % max_len)
    print("Average_patristic_distance\t%f" % (avg_distance_leaves/(nb_of_distances*scale)))

    branch_len_out.close()
    tree_info.close()

def ray_surveyor_config():

    print("TODO")


# Main #
if __name__ == "__main__":

    usage ='''
    python pyvolve-genseq.py <tree.nwk> <seq-size> [<scale> default=1 (no scale)]
    '''
    if len(sys.argv) < 3:
        sys.exit(usage)

    tree_f = sys.argv[1]
    outfiles = tree_f
    size = sys.argv[2]
    scale = 1
    scale = float(sys.argv[3]) if len(sys.argv) > 3 else None

    print("Reading tree..")
    my_tree = pyvolve.read_tree(file = tree_f, scale_tree=scale)
    my_model = pyvolve.Model("nucleotide")
    my_partition = pyvolve.Partition(models = my_model, size = int(size))

    print("Simulating sequences..")
    my_evolver = pyvolve.Evolver(tree = my_tree, partitions = my_partition)
    my_evolver(ratefile = "%s.%s.ratefile.txt" % (outfiles, size),
               infofile = "%s.%s.infofile.txt" % (outfiles, size),
               seqfile = "%s.%s.seqfile.fasta" % (outfiles, size) )

    print("Tree info..")
    tree_distances_info(tree_f, scale, int(size))

    print("Running ANI on sequences..")
    pyani_seq("%s.%s.seqfile.fasta" % (outfiles, size))


