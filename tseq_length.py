#! /usr/bin/env python

#This script just makes a list of the lengths of the fasta sequences in a file

import sys #We want to be able to get information from the command line.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels

Usage = '''
tseq_length.py version 1.0 9 May 2014 Abby Moore
This script calculates the length of sequences in fasta files.  They can be aligned 
i.e., with gaps, or not.
tseq_length.py [name of the file with the list of alignments] [folder containing the alignments]
[outfile]
'''

if len(sys.argv) < 4:
	print Usage
else:
	ListFile = sys.argv[1]
	AlignFolder = sys.argv[2]
	OutFileName= sys.argv[3]
	
FileList = [ ] #The list of files from ListFile
LengthDict = defaultdict(dict) #The dictionary of sequence lengths

InFile = open(ListFile, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	FileList.append(Line)
InFile.close()

print("The names of %d alignments were read from the file %s." % (len(FileList), ListFile))

for FileName in FileList:
	InFileName = AlignFolder + FileName
	SeqName = ""
	Seq = ""
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		if Line[0] == ">":
			if SeqName != "":
				SeqLen = 0
				for Char in Seq:
					if Char != "-":
						SeqLen += 1
				LengthDict[FileName][SeqName] = SeqLen
				Seq = ""
			SeqName = Line[1:]
		else:
			Seq += Line
	#Now we need to do the last sequence
	SeqLen = 0
	for Char in Seq:
		if Char != "-":
			SeqLen += 1
	LengthDict[FileName][SeqName] = SeqLen
	InFile.close()

print ("The length information from the sequences in %d files was read from the folder %s." % (len(LengthDict), AlignFolder))

OutList = [ ]
for FileName in FileList:
	SeqList = LengthDict[FileName].keys()
	SeqList.sort()
	for SeqName in SeqList:
		Line = FileName+"\t"+SeqName+"\t"+str(LengthDict[FileName][SeqName])+"\n"
		OutList.append(Line)

OutFile = open (OutFileName, 'w')
for Line in OutList:
	OutFile.write(Line)
OutFile.close()

print ("The sequence lengths were written to the file %s." % (OutFileName))
