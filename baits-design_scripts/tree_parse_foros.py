#! /usr/bin/env python

#tree_parse.py version 1.0 4 April 2014 Abby Moore
#This script is supposed to take the output from tree_parse.r, make the new alignments, 
#and save them in the appropriate folders.

#The script itself is meant to run on my laptop, but the muscle script it outputs
#is formatted for Oscar.  This script is quite fast and can go through the entire set
#of groups (originally 5022, but then split, so more than that now) in less than 30 minutes.
#This makes it easier than dealing with writing all of the different sbatch lines.

#Infile format:
#Note that first line is some kind of header.
'''[0] filename:  /home/abby/transcriptomes/newgroups/nGroup_10000a.fa
[1] Group Number: 10000
[2] can be "split" (if it was further split) or "ready" (if it is ready for final analysis)
[2] fate of original tree: split
[3] and [4] can be "continue" or "rejected, ...." (with information on why it was rejected)
[3] (if applicable), fate of first half of tree: rejected, NumCact: 10 NumMoll: 1 NumPort: 0
[4] (if applicable), fate of second half of tree: rejected, NumCact: 4 NumMoll: 0 NumPort: 1
'''
#First group saved in GroupNum_Alist.txt; second group saved in GroupNum_Blist.txt.

from collections import defaultdict #We will want to make dictionaries with multiple levels.
import re #We will want to write regular expressions
import sys #We want to be able to output error messages to the screen.
import os #We need to be able to talk to the operating system

if len(sys.argv) <4:
	print Usage
else:
	DirectoryIn = sys.argv[1] #where we find the lists of taxa and sequence files
	DirectoryOut = sys.argv[2] #where we put the new sequence files
	DirectoryReady = sys.argv[3] #where we put the finished sequence files that are ready for final analysis
	InFileName = sys.argv[4] #the file output by R

InFileCombined = DirectoryIn + InFileName

#The dictionaries and lists to be filled with information about the sequences and output
#files.
FateDict = { }
RList = [ ]
MuscleList = [ ]

#Reading the output from R, so we know what to do with the various files.
LineNum = 0
Finished = 0
AnalyzeFurther = 0
StopAnalysis = 0
InFile = open(InFileCombined, 'rU')
for Line in InFile:
	if LineNum > 0:
		Line = Line.strip('\r').strip('\n').split('\t')
		GroupNum = Line[1]
		#If a group of sequences is considered to be ready for the final tree building
		#(and does not need to be further split), then its alignment is moved to the 
		#directory of finished alignments.
		if Line[2] == "ready":
			OutLine = "cp "+DirectoryIn+"nGroup_"+GroupNum+"a.fa "+DirectoryReady+"nGroup_"+GroupNum+"r.fa"
			os.popen(OutLine,'r')
			FateDict[GroupNum] = "ready"
			Finished += 1
		#If not, we keep track of what to do with the two subtrees.
		elif Line[2] == "split":
			if Line[3] == "continue":
				FateDict[GroupNum+ "_A"] = "continue"
				AnalyzeFurther += 1
			else:
				FateDict[GroupNum+"_A"] = "stop"
				StopAnalysis += 1
			if Line[4] == "continue":
				FateDict[GroupNum+"_B"] = "continue"
				AnalyzeFurther += 1
			else:
				FateDict[GroupNum+"_B"] = "stop"
				StopAnalysis += 1
		else:
			print("Error with Group %s!!\n" % (GroupNum))
	LineNum += 1
InFile.close()

print ("%d lines were read from the input file %s.\n" % (LineNum, DirectoryIn+InFileName))
print ("%d sets of sequences are ready for probe design.  %d need further analysis.\n" % (Finished, AnalyzeFurther))
print ("%d were rejected.\n" % (StopAnalysis))

#making new sequence files for the groups of sequences that will be analyzed further.

for GroupNum in FateDict.keys():
	if FateDict[GroupNum] == "continue":
		RootGroupNum = GroupNum[:-2]
		#the file with all sequences (for both subgroups)
		SeqFile = DirectoryIn+"nGroup_"+RootGroupNum+".fa"
		#the file with the names of the sequences for that subgroup
		NamesFile = DirectoryIn+"nGroup_"+GroupNum+"list.txt"
		TempSeqDict = { }
		GroupSeqList = [ ]
		OutFileName = DirectoryOut+"nGroup_"+GroupNum+".fa"
		#The set of all of the sequences (from both subtrees) is read.
		InFile = open(SeqFile, 'rU')
		for Line in InFile:
			Line = Line.strip('\r').strip('\n')
			if Line[0] ==  ">":
				SeqName = Line[1:]
			else:
				Seq = Line
				TempSeqDict[SeqName] = Seq
		InFile.close()
		#Then the file with the names of the sequences that belong in the subtree
		#of interest is read, and those sequences are added to a separate list.
		InFile = open(NamesFile, 'rU')
		for Line in InFile:
			SeqName = Line.strip('\r').strip('\n')
			Seq = TempSeqDict[SeqName]
			OutLine = ">"+SeqName+"\n"+Seq+"\n"
			GroupSeqList.append(OutLine)
		InFile.close()
		#This is written to a file for alignment.
		OutFile = open(OutFileName, 'w')
		for Line in GroupSeqList:
			OutFile.write(Line)
		OutFile.close()
		#print ("For group %s, %d sequences were written to the outfile %s.\n" % (GroupNum, len(GroupSeqList), OutFileName))
		#A line is added to the file for R, so it knows to analyze the alignment.
		Line = DirectoryOut+"nGroup_"+GroupNum+"a.fa\t"+GroupNum+"\n"
		#alternative line, if I am running R on Oscar:
		#Line = "nGroup_"+GroupNum+"a.fa\t"+GroupNum+"\n"
		RList.append(Line)
		#And a line is added to the file for muscle, so it knows to make the alignment.
		Line = "muscle -in nGroup_"+GroupNum+".fa -out nGroup_"+GroupNum+"a.fa -maxmb 2000 -maxhours 2\n"
		MuscleList.append(Line)

#writing files for R and for muscle
OutFileName = DirectoryOut+"RList.txt"
OutFile = open(OutFileName, 'w')
for Line in RList:
	OutFile.write(Line)
OutFile.close()

#The muscle script is more complicated, because it is a script to run on Oscar.
OutFileName = DirectoryOut+"musclescript.sh"
OutFile = open(OutFileName, 'w')
Line = "#! /bin/bash\n#SBATCH -J Muscle\n#SBATCH -t 20:00:00\n#SBATCH -n 1\n\n"
OutFile.write(Line)
Line = "module load muscle\n"
OutFile.write(Line)
for Line in MuscleList:
	OutFile.write(Line)
OutFile.close()
			



