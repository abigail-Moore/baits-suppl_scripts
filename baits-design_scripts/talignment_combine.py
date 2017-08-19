#! /usr/bin/env python

#talignment_combine.py version 1.0 12 May 2014 Abby Moore
#This script reads many alignments in separate files and combines them into one file,
#while labeling the sequences with the alignment they were in, so that other sequences
#can be BLASTed against the alignment.

import sys #We want to be able to get data from the command line
import re #regular expressions

Usage = '''
talignment_combine.py [file that lists the alignments of interest] [alignmentfolder]
[name for the outfile]
The file listing the alignments of interest and the outfile name can be either the entire
path or just the file names if you are in that directory.  However, the file names in the 
file must only be the filename, not the full path.
It expects the sequence files to be in fasta format with no header lines.
'''

if len(sys.argv) < 4:
	print Usage
else:
	ListFile = sys.argv[1]
	SeqDirect = sys.argv[2]
	OutFileName= sys.argv[3]

FileList = [ ] #The list of files
SeqsOut = [ ] #The sequences to be written to the output file
NameRe = r"(\w+)\..+"
NameSub = r"\1"

#reading the list of files
InFile = open(ListFile, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	FileList.append(Line)
InFile.close()

#reading the sequences from each file
for FileName in FileList:
	InFile = open(SeqDirect+FileName, 'rU')
	FilePart = re.sub(NameRe, NameSub, FileName)
	SeqName = ""
	for Line in InFile:
		Line = Line.strip('\r').strip('\n')
		if Line[0] == ">":
			#first add the previous sequence (if there is one)
			if SeqName != "":
				OutLine = ">"+SeqName+"\n"+Seq+"\n"
				SeqsOut.append(OutLine)
			#then get a new sequence name
			SeqName = FilePart+"+"+Line[1:]
			#and clear the previous sequence
			Seq = ""
		else:
			Seq += Line
	#adding the last sequence
	OutLine = ">"+SeqName+"\n"+Seq+"\n"
	SeqsOut.append(OutLine)
	InFile.close()

#writing the sequences to a file
OutFile = open(OutFileName, 'w')
for Line in SeqsOut:
	OutFile.write(Line)
OutFile.close()
