'''
Created on Oct 12, 2020
Modified Nov 6, 2020

Flow:
1. List genes from gtf file (start, stop, strand, gene id)
2. List SNPs from input file. If SNP crosses thershold, mark it has a critical hit.
3. Count number of SNPs in each gene. Also calculate number of critical hits per gene and ratio of critical to all hits. Also calculate normalized score based on gene length.
4. Output .csv file with gene id, no. of hits, no. of critical hits, relative critical hits (cricital hits divided by total hits), hits normalized by length, gene start, gene stop and gene length

Requires biopython (python 3)
See --help for arguments

@author: Ville Hoikkala

'''
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

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("SNPmapper_log.txt", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

sys.stdout = Logger()

class SNP:
    position = 0
    freq = 0
    passed = 0

class ORF:
    def __init__(self):
        self.start = 0
        self.stop = 0
        self.ID = ""
        self.strand = ""
        self.length= 0
        self.hits = 0
        self.critical_hits = 0
        self.relative_critical_hits = 0
        self.normalizedHitsLength = 0

def getSNPpositions(file,thresholdMode,thresholdValue):
    parsed = list(csv.reader(open(file, 'r'), delimiter='\t')) #opens file as tab-delimited
    #spacerList = [spacer() for i in range(numberOfEntries)] #Create list to store spacers as objects
    SNPList = [] #Create list to store spacers as objects
    SNPCount = 0 #keep count of snips
    if threshold == "no":
        for line in parsed:
            SNPList.append(int(line[1]))
            SNPCount = SNPCount + 1
    if threshold == "yes":
        for line in parsed:
            #if line[]
            SNPList.append(int(line[1]))
            SNPCount = SNPCount + 1
    print("    Number of SNPs: " + str(SNPCount))
    return SNPList

def getORFs(csvfile,mode,contig):
    file = csvfile
    rowCounter = 0
    ORFcount = 0
    ORFlist = []
    for row in file: #Create an ORF object from each row
        if row[2] == mode and row[0] == contig:
            id = row[8]
            start = row[3]
            stop = row[4]
            length = int(stop)-int(start)
            direction = row[6]
            ORFlist.insert(ORFcount, ORF())
            #print("ID: " + id + ", start: " + start+ ", stop: " + stop + ", direction: " +direction)
            setattr(ORFlist[ORFcount], "id", id)
            setattr(ORFlist[ORFcount], "start", int(start))
            setattr(ORFlist[ORFcount], "stop", int(stop))
            setattr(ORFlist[ORFcount], "strand", str(direction))
            setattr(ORFlist[ORFcount], "length", int(length))
            ORFcount = ORFcount + 1
        rowCounter = rowCounter + 1
    print ("    Number of ORFs: " + str(len(ORFlist)))
    return ORFlist

def getSNPs(file,threshold,snp_freq_column,threshold_direction,contig):
    parsed = list(csv.reader(open(file, 'r'), delimiter='\t')) #opens file as tab-delimited
    SNPList = [] #Create list to store spacers as objects
    SNPCount = 0 #keep count of snips
    rowCounter = 0
    for row in parsed:
        if row[0] == contig:
            position = row[1]
            freq = row[int(snp_freq_column)]
            if threshold_direction == "under":
                #print("freq:" + str(freq))
                #print("threshold:" + str(threshold))
                if freq < threshold:
                    #print("BING" + str(freq))
                    passed = 1
                else:
                    passed = 0
            if threshold_direction == "over":
                if freq > threshold:
                    passed = 1
                else:
                    passed = 0
            SNPList.insert(SNPCount, SNP())
            setattr(SNPList[SNPCount], "position", int(position))
            setattr(SNPList[SNPCount], "freq", float(freq))
            setattr(SNPList[SNPCount], "passed", int(passed))
            SNPCount = SNPCount + 1
            rowCounter = rowCounter + 1
    if len(SNPList) == 0: print("No SNPs in given contig. Are you sure you gave the same contig for genes and SNPs?")
    return SNPList




def filterByContig(allRows,contig):
    filtered_ORFs = []
    for row in allRows:
        if row[0] == contig:
            filtered_ORFs.append(row)
    return filtered_ORFs

#python SNP_mapper_2.py --genes dada/monCan3F9.scf.v1.braker.chrom5.gtf --snps dada/hyb.txt  --output dada/out.test --output out.txt --mode gene --threshold_direction under --contig monCan3F9.scf.00004
# --snp_freq 2 --threshold 0.05
""" *** PROGRAM STARTS HERE *** """
print("")
print("")
print("-------------------")
print("Starting SNP mapper v 2.0")
print("Date " + str(datetime.datetime.now()))
print("")

print("INPUT:")

parser = argparse.ArgumentParser()
parser.add_argument("-g","--genes", help="input gtf file with predicted gene regions")
parser.add_argument("-s","--snps", help="file containing SNPs")
parser.add_argument("-o","--output", help="output filename")
parser.add_argument("-m","--mode", help="either gene or transcript")
parser.add_argument("-d","--threshold_direction", help="under or over given --threshold")
parser.add_argument("-c","--contig", help="contig name (for example monCan3F9.scf.00004)")
parser.add_argument("-sn","--snp_freq_column", help="column number for SNP frequency (first column is 0)")
parser.add_argument("-t","--threshold", help="frequency limit, e.g. 0.05")

args = parser.parse_args()

ORFfile = args.genes
SNPfile = args.snps
outputFilename = args.output
mode = args.mode
threshold_direction = args.threshold_direction
threshold = args.threshold
contig = args.contig
snp_freq_column = args.snp_freq_column

upOutfileO = outputFilename+"_ORF_stats.csv" #actual output filename for ORF stats

print("Gene input file: " + ORFfile)
print("SNP input file: " + SNPfile)
print("Output file: " + upOutfileO)
print("Mode: " + mode)
print("Threshold: " + threshold)
print("Threshold direction: " + threshold_direction)
print("")
print("RESULTS:")

#Open file
with open(ORFfile, 'rt') as csvfile:
    ORFreader = csv.reader(csvfile, delimiter = '\t')
    ORFs = list(ORFreader)

#Filter genes by contig
# print("Genes before filtering by contig: " + str(len(ORFs)))
# filtered_ORFs = filterByContig(ORFs,contig)
# print("Genes after filtering by contig:" + str(len(filtered_ORFs)))

#Filter SNPs by contig
# print("SNPs before filtering by contig: " + str(len(ORFs)))
# filtered_O = filterByContig(ORFs,contig)
# print("SNPs after filtering by contig:" + str(len(filtered_ORFs)))

#Make a list with all ORF objects (only from desired contig)
ORFobjects = getORFs(ORFs,mode,contig)

#Get SNPs getSNPs(SNPs,threshold,snp_freq_column,threshold_direction)
SNPobjects = getSNPs(SNPfile,threshold,snp_freq_column,threshold_direction,contig)
print(len(SNPobjects))
SNPs_orig = SNPobjects.copy()
SNPs_on_ORFs = 0

# #Get SNPs passing threshold
# SNPs_passed = getSNPpositions(SNPfile,"yes",threshold)
# SNPs_orig_passed = SNPs.copy()
# SNPs_on_ORFs_passed = []

progress = 0
intervalCounter = 0
intervalMemory = 0
interval = 0.01
totalORFs = len(ORFobjects)
for ORF in ORFobjects:
    for snp in SNPobjects:
        if snp.position >= ORF.start and snp.position <= ORF.stop:
            ORF.hits = ORF.hits + 1
            #print("SNP pass value:" + str(snp.passed))
            ORF.critical_hits = ORF.critical_hits + snp.passed
            SNPobjects.remove(snp)
            SNPs_on_ORFs = SNPs_on_ORFs + 1
    progress = progress +1
    intervalCounter = intervalCounter + 1
    if intervalCounter/totalORFs >= interval:
        intervalMemory = intervalMemory + interval
        print("\r" + str(round(Decimal(intervalMemory*100),0)) + "%", sep=' ', end='', flush=True)
        intervalCounter = 0


# #Calculate normalized hit score based on gene length and the proportion of hits crossing threshold
for o in ORFobjects:
    norm_score = round(Decimal(o.hits/o.length*10000),2)
    critical_hits =  o.critical_hits #hits passing threshold
    if o.hits > 0:
        relative_critical_hits =  round(Decimal(o.critical_hits/o.hits),2)
        setattr(o, "relative_critical_hits", relative_critical_hits)
    setattr(o, "norm_score", norm_score)
    setattr(o, "critical_hits", critical_hits)

# #Sort the ORF list
ORFobjects.sort(key=lambda x: x.relative_critical_hits, reverse=True)

# #Output csv file with the following header row
outList = []
outTable = [["gene_id","hits","critical_hits","relative_critical_hits","norm_score_length","gene_start","gene_stop","length"]]

# #Print each ORF into the csv file
for o in ORFobjects:
    values = [o.id,str(o.hits),str(o.critical_hits),str(o.relative_critical_hits),str(o.norm_score),str(o.start),str(o.stop),str(o.length)]
    outTable.append(values)
#
with open(upOutfileO, "w") as outfile:
    writer = csv.writer(outfile)
    writer.writerows(outTable)

#Stats for console
totalSNPs = len(SNPs_orig)
print("    Number of SNPs on genes: " + str(SNPs_on_ORFs) + " (" + str(round(Decimal((SNPs_on_ORFs/totalSNPs)*100),2)) + "%)")
print("    Number of intergenic SNPs : " + str(len(SNPobjects)) + " (" + str(round(Decimal((len(SNPobjects)/totalSNPs)*100),2)) + "%)")
print("")
print ("DONE!")
print ("----------------")
print("")
print("")
