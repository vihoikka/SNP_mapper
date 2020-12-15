from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import csv
from decimal import Decimal
from math import factorial
import argparse
import sys
import datetime

parser = argparse.ArgumentParser()
parser.add_argument("-g","--genes", help="input gtf file with predicted gene regions", required = True)
parser.add_argument("-l","--locations", help="locations to search genes from", required = True)
parser.add_argument("-o","--outfile", help="output filename", required = True)
args = parser.parse_args()


geneFile = args.genes
locFile = args.locations
outputfile = args.outfile

with open(geneFile, 'r') as csvfile:
    ORFreader = csv.reader(csvfile, delimiter = '\t')
    ORFs = list(ORFreader)

with open(locFile, 'r') as csvfile:
    locReader = csv.reader(csvfile, delimiter = '\t')
    locations = list(locReader)

    #parsed = list(csv.reader(open(file, 'r'), delimiter='\t')) #opens file as tab-delimited


chosenGenes = []
genecounter = 0

for loc in locations:
    loc_start = int(loc[1])
    loc_stop = int(loc[2])
    #print(loc_start)
    #print(loc_stop)
    for g in ORFs:
        gene_start = int(g[3])
        gene_stop = int(g[4])
        #print("gene start: " + str(gene_start) + ". Stop: " + str(gene_stop))
        if (((gene_start < loc_start and gene_stop > loc_start) or   (gene_start > loc_start and gene_stop > loc_start and gene_stop < loc_stop) or  (gene_start < loc_stop and gene_stop > loc_stop))  and g[0] == loc[0]):
            #if g not in chosenGenes:
            chosenGenes.append(g)
            #print("added " + str(g))
            #ORFs.remove(g)
            genecounter =+ 1

with open(outputfile, "w") as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerows(chosenGenes)
