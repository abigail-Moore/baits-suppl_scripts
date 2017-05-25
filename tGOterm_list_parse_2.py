#! /usr/bin/env python

#This is supposed to find the GO terms of interest and print the lines of the GOspread.txt 
#file that have those GO terms, along with the meaning of the GO term.

import sys #We want to be able to get file names from the command line
from collections import defaultdict #We want to be able to make dictionaries with multiple levels.

Usage = '''
tGOterm_list_parse_2.py 7 May 2014 Abby Moore
This script expects the following input:
tGOterm_list_parse.py [file with GO terms of interest] [GOspread.txt file] 
[file with information about the alignments] [file with the AT locus names]
[prefix for the outfiles]
The files can either be the complete path to those files, if you are not in the directory
that contains those files or just the file names if you are in that directory.
'''

if len(sys.argv) < 6:
	print Usage
else:
	InFileGO = sys.argv[1]
	InFileGOspread = sys.argv[2]
	InFileALinfo = sys.argv[3]
	InFileATinfo = sys.argv[4]
	OutFilePre = sys.argv[5]

#The expected format of the file of GO terms (tab-delimited):
#GO:0009769 [0]
#photosynthesis, light harvesting in photosystem II [1]

#The expected format of the GOspread.txt file (tab-delimited):
#0000001 [0]: GO term
#9350_B_A_A [1]: alignment name
#Univ [2]: code (Univ for in every AT gene of that orthogroup, All for in every AT gene
#of that alignment, Some for in at least one of the AT genes of that alignment)
#0051276,0005739,0003876,..... [3]: other GO terms for that alignment

#the expected format of the ALinfo file (tab-delimited):
#10011_B [0]: alignment name
#CUVY_2099504:186 [1]: Mollugo sequences and lengths
#Pe-bleo_69624.0_3:2477,Pe-bleo_69624.0_2:2543..... [2]: Cactaceae sequences and lengths
#An-fil_38684.0_4:1977,An-fil_38684.0_2:2137..... [3]: other Portulacineae sequences and lengths
#AT1G51965 [4]: Arabidopsis sequences
#0005737,0034641... [5]: universal GO terms for that orthogroup
# (if present) [6]: GO terms present in all of the AT sequences in the alignment
# (if present) [7]: GO terms present in some of the AT sequences in the alignment

#the expecect format of the ATinfo file (tab-delimited):
#AT1G01010 [0]: AT locus name
#NAC001 [1]: abbreviated gene name
#NAC domain containing protein 1 [2]: long gene name

#The outfile will have the first three columns [0-2] of GOspread, the definition of the GO
#term (column 1 from GO terms of interest), and the last column of GO spread.

#the outfiles:
OutFileGO_AL = OutFilePre+"GO_to_AL.txt"
OutFileAL_GO = OutFilePre+"AL_to_GO.txt"
OutFileGO = OutFilePre+"GO.txt"

#the dictionaries and things we fill out
GOInterest = { } #from InFileGO
GOSpread = defaultdict(dict) #from InFileGOspread

ALGODict = { }
ATNameDict = { } #from InFileATinfo
ALATDict = { } #from InFileALinfo

InFile = open(InFileGO, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	GOTerm = Line[0][3:]
	GOInterest[GOTerm] = Line[1]
InFile.close()

GOList = GOInterest.keys()
GOList.sort()

print ("%d GO terms of interest were read from infile %s." % (len(GOInterest), InFileGO))

InFile = open(InFileGOspread, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	GOTerm = Line[0]
	ALName = Line[1]
	try:
		GOSpread[GOTerm][ALName]['Code'] = Line[2]
	except KeyError:
		GOSpread[GOTerm][ALName] = defaultdict(dict)
		GOSpread[GOTerm][ALName]['Code'] = Line[2]
	#I am being lazy here and not splitting the list, so it is just a single string.
	GOSpread[GOTerm][ALName]['GOList'] = Line[3]
InFile.close()

print ("Information about %d GO terms was read from infile %s."% (len(GOSpread.keys()), InFileGOspread))

InFile = open(InFileALinfo, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	ALName = Line[0]
	ATList = Line[4].split(',')
	ALATDict[ALName] = ATList
InFile.close()

print("Information about %d alignments was read from infile %s." % (len(ALATDict), InFileALinfo))

InFile = open(InFileATinfo, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	ATName = Line[0]
	GeneName = Line[2]
	ATNameDict[ATName] = GeneName
InFile.close()

print("The names of %d Arabidopsis loci was read from infile %s." % (len(ATNameDict), InFileATinfo))

for GOTerm in GOList:
	for ALName in GOSpread[GOTerm].keys():
		try:
			ALGODict[ALName].append(GOTerm)
		except KeyError:
			ALGODict[ALName] = [GOTerm]

OutList = [ ]
for GOTerm in GOList:
	for ALName in GOSpread[GOTerm].keys():
		Line = GOTerm+"\t"+ALName+"\t"+GOSpread[GOTerm][ALName]['Code']+"\t"
		Line += GOInterest[GOTerm]+"\t"+GOSpread[GOTerm][ALName]['GOList']+"\n"
		OutList.append(Line)

OutFile = open(OutFileGO_AL, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()

print("Information about the alignments in the GO terms of interest was written to the file %s."% (OutFileGO_AL))

ALList = ALGODict.keys()
ALList.sort()

OutList = [ ]
for ALName in ALList:
	Line = ALName+"\t"
	for GOTerm in ALGODict[ALName]:
		Line += GOTerm+": "+GOInterest[GOTerm]+" "+GOSpread[GOTerm][ALName]['Code']+"\t"
	for ATName in ALATDict[ALName]:
		try:
			Line += ATName+": "+ATNameDict[ATName]+"\t"
		except KeyError:
			Line += ATName+": no name\t"
	Line = Line[:-1]+"\n"
	OutList.append(Line)
OutFile = open(OutFileAL_GO, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()

print("Information about the alignments that have those GO terms was written to the file %s." % (OutFileAL_GO))

OutList = [ ]
for GOTerm in GOList:
	Line = GOTerm+"\t"+GOInterest[GOTerm]+"\t"+str(len(GOSpread[GOTerm].keys()))+"\t"
	if len(GOSpread[GOTerm].keys()) != 0:
		for ALName in GOSpread[GOTerm].keys():
			Line += ALName+"("+str(len(ALGODict[ALName]))+"), "
			Line = Line[:-1]
	Line = Line[:-1]+"\n"
	OutList.append(Line)
OutFile = open(OutFileGO, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()

print("Condensed information about the GO terms was written to the file %s." % (OutFileGO))
