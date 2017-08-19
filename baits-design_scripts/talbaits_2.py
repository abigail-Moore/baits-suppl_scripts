#! /usr/bin/env python

#talbaits.py version 2.0 3 June 2014 Abby Moore

#This script takes a list of sequences in alignments (as produced by talsplit.py or tseq_length.py),
#trims them accordingly, and makes a new file for each alignment that only has the sequences
#we want to be made into baits.  Then it writes a script to align the sequences.

#The expected format of the sequence list file (tab delimited):
'''
mdh5.fa [0]: alignment file name
scaffold-HURS-2047554-Mollugo_pentaphylla-mature_leaf [1]: sequence name
1611 [2]: sequence length
1 [350:1350] [3]: number of copies of the sequence desired (1 or 2), possibly followed by
	a range in parentheses if we do not want the whole sequence; if we do not want the sequence,
	this will not exist and [3] will be the comments
similar to BZMI, randomly chose [4]: comments
'''

import sys #We want to be able to get information from the command line.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels

Usage = '''
talbaits_2.py
This is a script to take a list of sequences in alignments, find those sequences, and write
them to alignments for which baits can be made.
talbaits_2.py [sequence list] [folder where alignments are found]
[folder where new alignments should be put]
It does not assume that the first and last files are/should be in the alignment folder.
'''

if len(sys.argv) < 4:
	print(Usage)
else:
	ListFileName = sys.argv[1]
	SeqFilePre = sys.argv[2]
	OutFilePre = sys.argv[3]


SeqListDict = defaultdict(dict) #This will list the sequences and their information grouped
	#by alignment
NumSeqs = 0
AlScript = [ ]

#adding a slash to the end of the file path name if it is not already there
if SeqFilePre[-1] != "/":
	SeqFilePre.append("/")
#and doing the same for the out file path
if OutFilePre[-1] != "/":
	OutFilePre.append("/")

#first to read the list of baits we want
InFile = open(ListFileName, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	#skipping the line if it is blank or a comment
	if Line != "" and Line[0] != "#":
		Line = Line.split('\t')
		#skipping the line if it is somehow too short or if we do not want the sequence
		if len(Line) > 3 and (Line[3][0] == "1" or Line[3][0] == "2"):
			AlName = Line[0]
			SeqName = Line[1]
			SeqLen = int(Line[2])
			#reading the range of the sequence we want, if we don't want the entire thing
			if len(Line[3]) > 1:
				RangeTemp = Line[3][1:]
				RangeTemp = RangeTemp.strip(" ").strip("[").strip("]").split(":")
				try:
					if RangeTemp[0] != "":
						CStart = int(RangeTemp[0])
					else:
						CStart = 0
				except IndexError:
					print Line
				try:
					if RangeTemp[1] != "":
						CEnd = int(RangeTemp[1])
					else:
						CEnd = SeqLen
				except IndexError:
					print Line
			else:
				CStart = 0
				CEnd = SeqLen
			NumSeqs += 1
			SeqListDict[AlName][SeqName] = defaultdict(dict)
			SeqListDict[AlName][SeqName]['SeqLen'] = SeqLen
			SeqListDict[AlName][SeqName]['CStart'] = CStart
			SeqListDict[AlName][SeqName]['CEnd'] = CEnd
InFile.close()

#print("InFile %s was read.  We want to design baits for %d loci.\n" % (ListFileName, len(SeqListDict)))
sys.stderr.write("InFile %s was read.  We want to design baits for %d loci.\n" % (ListFileName, len(SeqListDict)))

#print("We will use %d template sequences.\n" % (NumSeqs))
sys.stderr.write("We will use %d template sequences.\n" % (NumSeqs))

#Writing the start of the alignment script
Line = "#! /bin/bash\n\n"
AlScript.append(Line)

AlNum = 0
#making the alignment files
for AlName in SeqListDict:
	ToWrite = [ ]
	#reading the sequences from the file
	InFileName = SeqFilePre+AlName
	InFile = open(InFileName, 'rU')
	SeqsDictTemp = { }
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		if Line[0] == ">":
			SeqName = Line[1:]
			SeqsDictTemp[SeqName] = ""
		else:
			SeqsDictTemp[SeqName] += Line
	InFile.close()
	for SeqName in SeqListDict[AlName]:
		SeqTemp = SeqsDictTemp[SeqName]
		#if the file is of aligned sequences, then the gaps need to be removed.
		if "-" in SeqTemp:
			NewSeq = ""
			for base in SeqTemp:
				if base != "-":
					NewSeq += base
			SeqTemp = NewSeq
		#reading in the remaining data from SeqListDict
		SeqLen = SeqListDict[AlName][SeqName]['SeqLen']
		CStart = SeqListDict[AlName][SeqName]['CStart']
		CEnd = SeqListDict[AlName][SeqName]['CEnd']
		if len(SeqTemp) != SeqLen:
			sys.stderr.write("ERROR, the length of sequence %s is %d, but should be %d!!\n" % \
				(SeqName, len(SeqTemp), SeqLen))
		SeqTemp = SeqTemp[CStart:CEnd]
		Line = ">"+SeqName+"\n"+SeqTemp+"\n"
		ToWrite.append(Line)
	OutFileName = OutFilePre + AlName
	OutFile = open(OutFileName, 'w')
	for Line in ToWrite:
		OutFile.write(Line)
	OutFile.close()
	Line = "muscle3.8.31_i86linux64 -in "+OutFilePre+AlName+" -out "+OutFilePre+AlName[:-3]+"a.fa\n"
	AlScript.append(Line)
	AlNum += 1

if AlNum != len(SeqListDict):
	sys.stderr.write("ERROR!  %d files of sequences were written, but %d files should have been \
	written!!\n" % (AlNum, len(SeqListDict)))

#print("%d groups of were written to the directory %s.\n" % (AlNum, OutFilePre))
sys.stderr.write("%d groups of sequences were written to the directory %s.\n" % (AlNum, OutFilePre))

OutFileName = OutFilePre+"AlScript.sh"
OutFile = open(OutFileName, 'w')
for Line in AlScript:
	OutFile.write(Line)
OutFile.close()

#print("The script to align these sequences is called %s.\n" % (OutFileName))
sys.stderr.write("The script to align these sequences is called %s.\n" % (OutFileName))
