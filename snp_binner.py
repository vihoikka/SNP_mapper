# python snp_binner.py --bins 100 --genome_length 29432383 --snps data_4.csv

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

""" *** PROGRAM STARTS HERE *** """
print("")
print("")
print("-------------------")
print("Starting SNP binner v 0.1")
#print("Date " + str(datetime.datetime.now()))
print("")

print("INPUT:")

parser = argparse.ArgumentParser()
parser.add_argument("-s","--snps", help="input gtf file with predicted gene regions")
parser.add_argument("-b","--bins", help="file containing SNPs")
parser.add_argument("-le","--genome_length", help="output filename")
parser.add_argument("-o","--output", help="output filename")
args = parser.parse_args()

snpsFile = args.snps
bins = int(args.bins)
gl = int(args.genome_length)

binsize = round(gl/bins,0)
print("bin size " + str(binsize))

binSNPs = [0] * bins
print(len(binSNPs))

binDic = {}

binCounter = 1
while binCounter <= bins:
    binDic[binCounter] = 0
    binCounter += 1


with open(snpsFile, 'rt') as csvfile:
    snpReader = csv.reader(csvfile, delimiter = '\t')
    snps = list(snpReader)

print(len(snps))

binCounter = 0
for bin, hits in binDic.items():
    #print(binsize)
    binstart = binCounter * binsize + 1
    binstop = binstart + binsize - 1
    #print("Bin start: " + str(binstart))
    #print("Bin stop: " + str(binstop))
    for row in snps:
        #print(row[1])
        if int(row[1]) > binstart and int(row[1]) < binstop:
            binDic[bin] += 1
    binCounter += 1

binCounter = 0
for bin, hits in binDic.items():
    binstart = binCounter * binsize + 1
    binstop = binstart + binsize - 1
    #print("Bin start: " + str(binstart))
    #print("Bin stop: " + str(binstop))
    print("Bin: " + str(bin) + ". SNPs: " + str(hits))
    #print(str(hits))
    binCounter += 1

#for i in snps:
    #print(i)

#with open(snps)

#for i in binSNPs
