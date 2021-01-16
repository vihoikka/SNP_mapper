'''
Created on Nov 30, 2017

@author: ville
'''
import csv
import sys

averagedValues = []
numberOfValuesToBeAveraged = 400

scaffoldList = []
scaffoldDict= {}
origFile = []

print ("Doing things")

with open(sys.argv[1], 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter="\t")
   
    for row in reader:
        origFile.append(row)
   
    csvfile.seek(0)
   
    for row in reader:
        scaf = row[0]
        chrom = row[6]
        if scaf not in scaffoldList:
            scaffoldList.append(scaf)
            scaffoldDict[scaf] = chrom
    print (len(scaffoldDict))
   
    csvfile.seek(1)
    rowCounter = 0
   
    #this dictionary stores the summed values of many rows
    columnTotals = {"scaffold": "-", "scaffold_pos": 0, "flavohybr_freq": 0, "montanahybr_freq": 0, "montana_freq": 0, "flavo_freq": 0, "chromosome": "-"}
   
    #this list will (in the end) contain averaged values for the number of rows desired
    averagedValues = []
   
    scaffoldCounter = 0
   
    for scaffold, chromosome in scaffoldDict.items():
        print (str(scaffoldCounter)) + "/" + (str(len(scaffoldDict)))
        #print("hello")
        #start by going through each row in the dataset
        for row in origFile: #this means that in each iteration of the for -loop, the variable 'row' refers to the row currently being processed.
            #check how many rows we have scanned so far. In the beginning it's zero, so we may proceed
            if rowCounter < numberOfValuesToBeAveraged and row[0] == scaffold:
               
                #These five lines add the values of the current row into the dictionary columnTotals, which contains the summed values of multiple rows
                #Note that python treats each row as a list, with different columns being different elements of the list
                columnTotals["scaffold"] = scaffold
                columnTotals["scaffold_pos"] = columnTotals["scaffold_pos"] + float(row[1])
                columnTotals["flavohybr_freq"] = columnTotals["flavohybr_freq"] + float(row[2])
                columnTotals["montanahybr_freq"] = columnTotals["montanahybr_freq"] + float(row[3])
                columnTotals["montana_freq"] = columnTotals["montana_freq"] + float(row[4])
                columnTotals["flavo_freq"] = columnTotals["flavo_freq"] + float(row[5]) #note that float is similar to int, but for numbers with decimal values
                columnTotals["chromosome"] = chromosome
               
                rowCounter = rowCounter + 1
           
            #next, check if we have already gone through 100 rows. If we have, then we must...
            if rowCounter == numberOfValuesToBeAveraged:
                columnTotals["scaffold_pos"] = columnTotals["scaffold_pos"]/numberOfValuesToBeAveraged              
                columnTotals["flavohybr_freq"] = columnTotals["flavohybr_freq"]/numberOfValuesToBeAveraged
                columnTotals["montanahybr_freq"] = columnTotals["montanahybr_freq"]/numberOfValuesToBeAveraged
                columnTotals["montana_freq"] = columnTotals["montana_freq"]/numberOfValuesToBeAveraged  
                columnTotals["flavo_freq"] = columnTotals["flavo_freq"]/numberOfValuesToBeAveraged

                #Now, instead of the total, the columnTotals dictionary contains the AVERAGE value of 100 rows for each parameter
               
                rowCounter = 0 #reset rowCounter. In the next iteration of the for loop, this will start again from zero. Once reaching 100 again, it will start over again from zero because of this row
               
                #Finally, we need to copy the values from the columnTotals dictionary to a list. Note that this list also contains strings "scaffold" and "chrom" for all rows in its corresponding positions (couldn't take average from these)
                averagedRow = [scaffold, str(columnTotals["scaffold_pos"]), str(columnTotals["flavohybr_freq"]), str(columnTotals["montanahybr_freq"]), str(columnTotals["montana_freq"]), str(columnTotals["flavo_freq"]),  chromosome]
               
                #This list will then be appended to the averagedValues list which in the end will become the file that is saved. In averagedValues, each list element is another list. List within a list. Like Inception.
                averagedValues.append(averagedRow)
               
                #we must also reset the dictionary for the upcoming 100 entries. The same dictionary is therefore used multiple times, but emptied between 100 batches of rows. This saves nature.
                columnTotals = {"scaffold": "-", "scaffold_pos": 0, "flavohybr_freq": 0, "montanahybr_freq": 0, "montana_freq": 0, "flavo_freq": 0, "chromosome": "-"}
        csvfile.seek(0)
        scaffoldCounter = scaffoldCounter + 1
           
#Once we have the averagedValues list, in which each element contains the average values for all columns from 100 rows, we must save it to a csv file using tab (\t) as a delimiter. If this doesn't work, you can also change it to space (" ") or comma (",")
with open("averagedFile400.txt", "wb") as csvoutputfile:
    writer = csv.writer(csvoutputfile, delimiter="\t")
    writer.writerows(averagedValues)
   
print ("Things are done yay")
