#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# date:    	2016-10-16
# version: 	0.01


import sys
import os

# Main #
if __name__ == "__main__":

    usage = '''
    usage : raysurveyor-conf.py <fasta dir> <kmer length>
    '''

    if len(sys.argv) < 3:
        print(usage)
        sys.exit()

    dir = sys.argv[1]
    if dir[-1] == "/":
        # print("replacing trailing")
        dir = dir[0:-2]

    klen = sys.argv[2]
    outdir = dir+".survey.k"+klen

    with open(outdir+".conf","w") as fout:
        fout.write("-k "+klen+"\n")
        fout.write("-o "+outdir+"\n")
        fout.write("-run-surveyor"+"\n")
        for f in os.listdir(dir):
            fout.write("-read-sample-assembly %s %s" %
                       (f.replace(".fasta",""),
                        dir+"/"+f+"\n"))

