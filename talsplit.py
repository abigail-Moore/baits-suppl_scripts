#! /usr/bin/env python

#talsplit.py version 1.0 Abby Moore 20 May 2014

#This script just splits fasta alignment files according to the group numbers given in a
#text file.

'''
It expects input in the following format:
ppt.fa [0]: alignment name
scaffold-UNSW-2078418-Mollugo_nudicaulis-mature_leaf [1]: sequence name
127 [2]: sequence length (or anything else; this column is ignored)
2 [3]: group number (or x for no group)
'''

import sys #We want to be able to get data from the command line
import re #regular expressions
from collections import defaultdict #We want dictionaries with multiple levels.

Usage = '''
This script splits fasta files into groups of sequences.
talsplit.py [directory in which alignments are saved] [file with information about the alignments]
The information file should have the alignment name in the first column, the sequence
name in the second column, and the new group in the fourth column (with x for delete sequence).
'''

if len(sys.argv)<3:
	print(Usage)
else:
	SeqDirect = sys.argv[1]
	InFileAl = sys.argv[2]

AlDict = defaultdict(dict)
SeqLenDict = { }

#reading the alignment information file and saving the information to a default dictionary
InFile = open(InFileAl,'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	#ignores blank lines and comments
	if Line != "" and Line[0] != "#":
		Line = Line.split('\t')
		try:
			if Line[3] != "x" and Line[3] != "X":
				AlName = Line[0][:-3]
				SeqName = Line[1]
				SeqLen = Line[2]
				SeqGroup = Line[3]
				SeqLenDict[SeqName] = SeqLen
				try:
					AlDict[AlName][SeqGroup].append(SeqName)
				except KeyError:
					AlDict[AlName][SeqGroup] = [SeqName]
		except IndexError:
			print Line
InFile.close()

print("Information about %d alignments was read from the file %s." % (len(AlDict), InFileAl))

#going through each alignment
GroupList = [ ]
SeqList = [ ]
for AlName in AlDict.keys():
	InFileName = SeqDirect+AlName+".fa"
	InFile = open(InFileName, 'rU')
	#reading the sequences to a dictionary
	TempSeqs = { }
	for Line in InFile:
		Line = Line.strip('\r').strip('\n')
		if Line[0] == ">":
			SeqName = Line[1:]
			TempSeqs[SeqName] = ""
		else:
			TempSeqs[SeqName] += Line
	InFile.close()
	#writing the sequences to files for the various groups
	for SeqGroup in AlDict[AlName].keys():
		TempList = AlDict[AlName][SeqGroup]
		TempList.sort()
		AlDict[AlName][SeqGroup] = TempList
		GroupName = AlName+SeqGroup
		OutFileName = SeqDirect+GroupName+".fa"
		OutFile = open(OutFileName, 'w')
		for SeqName in AlDict[AlName][SeqGroup]:
			Line = ">"+SeqName+"\n"+TempSeqs[SeqName]+"\n"
			OutFile.write(Line)
			#and add information about that sequence to the list of sequences
			Line = GroupName+".fa\t"+SeqName+"\t"+SeqLenDict[SeqName]+"\n"
			SeqList.append(Line)
		OutFile.close()
		#and making a list of the groups
		GroupList.append(GroupName)

print("These alignments were split into a total of %d groups and saved to files with names of the form %s."\
	% (len(GroupList), OutFileName))

#Finally, writing a shell script to align and build trees for each of the groups:
ShellScript = ["#! /bin/bash\n"]
for GroupName in GroupList:
	Line = "muscle3.8.31_i86linux64 -in "+GroupName+".fa -out "+GroupName+"a.fa\n\
	MFAtoPHY.pl "+GroupName+"a.fa\nraxmlHPC -s "+GroupName+"a.fa.phy -n "+GroupName+" -m GTRCAT -p 1234\n"
	ShellScript.append(Line)
OutFileName = SeqDirect+"AlScript.sh"
OutFile = open(OutFileName, 'w')
for Line in ShellScript:
	OutFile.write(Line)
OutFile.close()

print("The script to analyze these new alignments is called %s." % (OutFileName)) 

#write the list of sequences to a file
OutFileName = SeqDirect+"sequencelist.txt"
OutFile = open(OutFileName, 'w')
for Line in SeqList:
	OutFile.write(Line)
OutFile.close()

print("Information about these sequences was written to the file %s." % (OutFileName))
		


