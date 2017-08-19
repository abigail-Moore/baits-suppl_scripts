#! /usr/bin/env python

#blast_dict.py version 1.0 11 March 2014 Abby Moore
#This script should take the transcriptomes of interest (Portulacineae, Beta vulgaris, and Mollugo)
#BLAST them against the groups of orthologous loci, and add the sequences that have good BLAST hits
#to the group containing the sequence they BLASTed to.

#The columns in the default BLASTN tabular output:
#Fields: query id [0], subject id [1], % identity [2], alignment length [3], mismatches [4], gap opens [5], q. start [6], q. end [7], 
#s. start [8], s. end [9], evalue [10], bit score [11]

import re #We will need the regular expressions module
import sys #We want to be able to output error messages to the screen and get input from the command line.
import gzip #We want to be able to open zipped files.
import os #We want to be able to talk to the shell.
from collections import defaultdict #We want to make dictionaries with multiple levels.

#input files:
#FileList = '/home/abby/transcriptomes/reference_transcriptomes/filelist.txt'
#DBFile = '/home/abby/transcriptomes/groups/Group_All.db'
#BOut = '/home/abby/transcriptomes/transcriptomes_cleaned/outfile.txt'
FileList = '/home/abby/transcriptomes/Mollugo/filelist.txt'
#DBFile = '/home/abby/transcriptomes/C4alignmentsPA/plusMoll/C4combal.db'
DBFile = '/home/abby/transcriptomes/C4alignmentsPA/plusMoll/mdhdb'
BOut = '/home/abby/transcriptomes/C4alignmentsPA/plusMoll/outfile.txt'


#output:
GroupDict = defaultdict(dict) #This will hold the sequences sorted by group number.
#OutFilePre = '/home/abby/transcriptomes/groups/Group_'
OutFilePre = '/home/abby/transcriptomes/C4alignmentsPA/plusMoll/'
OutFilePost = '.fa'
#OutFilePost = 'temp.fa'

#Making the list of input files:
InFileList = [ ]
NumFiles = 0
InFile = open(FileList, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	InFileList.append(Line)
	NumFiles += 1
InFile.close()
print ("The names of %d files were read from the file %s.\n" % (NumFiles, FileList))
sys.stderr.write("The names of %d files were read from the file %s.\n" % (NumFiles, FileList))

for InFileName in InFileList:
	print ("Analyzing the transcriptome %s.\n" % (InFileName))
	sys.stderr.write("Analyzing the transcriptome %s.\n" % (InFileName))
	#setting or resetting the dictionaries for this transcriptome
	TempGDict = { }
	TempSDict = { }
	OldGene = ""
	Seq = ""
	SeqNum = 0
	GeneNum = 0
	GroupedGene = 0
	#Write the BLAST command and send it to the operating system to run.
	Line = "blastn -db "+DBFile +" -query "+InFileName+" -out "+BOut+" -outfmt 6"
	os.popen(Line)
	#Open the first transcriptome:
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n')
		if Line[0] == ">":
			if Seq != "":
				TempSDict[SeqName] = Seq
				SeqNum += 1
			SeqName = Line.strip('>').split(' ')[0]
			GeneName = SeqName.split('.')[0]
			Seq = ""
		else:
			Seq += Line
	TempSDict[SeqName] = Seq
	print ("Sequences of %d genes were read from this transcriptome.\n" % (SeqNum))
	sys.stderr.write("Sequences of %d genes were read from this transcriptome.\n" % (SeqNum))
	BlastRes = open(BOut, 'rU')
	for Line in BlastRes:
		Line = Line.strip('\n').strip('\r').split('\t')
		Gene = Line[0]
		if Gene != OldGene:
			Group = Line[1].split("_")[0]
			TempGDict[Gene] = Group
			OldGene = Gene
			Score = float(Line[11])
			GeneNum += 1
		#I wanted this next part is as a check, but it if it does not appear to ever happen after more trials, it should
		#probably be commented out to speed things up.
		elif Score < float(Line[11]):
			sys.stderr.write("Error, the BLAST hits for gene %s are out of order.\n" % (Gene))
			sys.stderr.write("old score: %f, new score: %s\n" % (Score, Line[11]))
	print ("%d of these genes had BLAST hits to the reference sequences.\n" % (GeneNum))
	sys.stderr.write("%d of these genes had BLAST hits to the reference sequences.\n" % (GeneNum))
	for Gene in TempGDict.keys():
		GroupDict[TempGDict[Gene]][Gene] = TempSDict[Gene]
		GroupedGene += 1
	print ("%d genes were added to the dictionary of grouped sequences.\n" % (GroupedGene))
	sys.stderr.write("%d genes were added to the dictionary of grouped sequences.\n" % (GroupedGene))
	InFile.close()

GroupNum = 0
SeqNum = 0
for Group in GroupDict.keys():
	SeqstoWrite = [ ]
	for Gene in GroupDict[Group].keys():
		Line = ">"+Gene+"\n"+GroupDict[Group][Gene]+"\n"
		SeqstoWrite.append(Line)
		SeqNum += 1
	OutFileName = OutFilePre+Group+OutFilePost
	OutFile = open(OutFileName, 'a')
	for Line in SeqstoWrite:
		OutFile.write(Line)
	OutFile.close()
	GroupNum += 1

print ("The sequences for %d sequences in %d groups were appended to the end of their respective group files.\n" % (SeqNum, GroupNum))
sys.stderr.write("The sequences for %d sequences in %d groups were appended to the end of their respective group files.\n" % (SeqNum, GroupNum))
