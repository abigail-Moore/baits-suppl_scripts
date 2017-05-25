#! /usr/bin/env python

#tgroup_filter.py version 2.0 31 March 2014 Abby Moore
#This script is supposed to look through the files of sequences that belong to each 
#group.  There should be two files per group.  One, labeled Group_xxx, includes the
#model organism sequences, while the other, labeled Group_xxxtemp.fa, includes the
#Portullugo sequences.
#First, it looks through the file with Portullugo sequences to determine if there is
#at least one Mollugo and one Portulacineae sequence.  If that is true, then it processes
#everything further, combining the two files into one and shortening the names.
#osc is modified to provide an output muscle file that Oscar can execute.

from collections import defaultdict #We will want to make dictionaries with multiple levels.
import re #We will want to write regular expressions
import sys #We want to be able to output error messages to the screen (maybe?).

#The files:
InFileList = '/home/abby/transcriptomes/groups/groupfilelist.txt'
#InFileList = '/home/abby/transcriptomes/groups/groupfilelisthead.txt'
PortNames = '/home/abby/transcriptomes/transtaxa2.txt'
PortCat = '/home/abby/transcriptomes/transtaxacat.txt'
InFilePre = '/home/abby/transcriptomes/groups/'
OutFilePre = '/home/abby/transcriptomes/newgroups/nGroup_'

#The dictionaries:
FileList = [ ]
GroupStats = defaultdict(list)
TaxaDict = { }
CatDict = { }
GroupFileList = [ ]
GroupFate = { }
NumFailed = 0
NumWorked = 0
ListAccepted = [ ]

#Regular expressions to find the group number from the file name:
ReFileName = r"Group_(\d+)temp\.fa"
SubGroupNum = r"\1"
SubFileName = r"Group_\1\.fa"

#Making a list of files that have Portullugo sequences.
InFile = open(InFileList, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	if Line.endswith('temp.fa'):
		FileList.append(Line)
InFile.close()

print ("Infile %s read.\n" % (InFileList))

#Making a dictionary to tell which transcriptomes belong to which taxa.
InFile = open(PortNames, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	TaxaDict[Line[1]] = Line[0]
InFile.close()

print ("Infile %s read.\n" % (PortNames))

#Making a dictionary to tell which transcriptomes belong to which groups.
InFile = open(PortCat, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	CatDict[Line[0]] = int(Line[1])
InFile.close()

print ("Infile %s read.\n" % (PortCat))

NumFiles = 0
NumFilesOut = 0
for File in FileList:
	GroupNum = re.sub(ReFileName,SubGroupNum,File)
	InFileName = InFilePre + File
	InFile = open(InFileName, 'rU')
	#resetting the various variables/lists
	SeqDictTemp = { }
	NameDictTemp = { }
	NumMoll = 0
	SeqLenMoll = [ ]
	NumBeta = 0
	SeqLenBeta = [ ]
	NumPort = 0
	SeqLenPort = [ ]
	NumCact = 0
	SeqLenCact = [ ]
	NumCary = 0
	SeqLenCary = [ ]
	#reading the sequences in
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		if Line[0] == ">":
			SeqName = Line[1:]
		else:
			Seq = Line
			SeqDictTemp[SeqName] = Seq
	InFile.close()
	#shortening the names
	for SeqName in SeqDictTemp.keys():
		SSeqName = SeqName.split("-")
		#The Mollugo names are very long and will be shortened to the
		#transcriptome code and transcript number.
		if SSeqName[0] == "scaffold":
			NewName = "M_"+SSeqName[1]+"_"+SSeqName[2]
			NameDictTemp[SeqName] = NewName
			NumMoll += 1
			SeqLenMoll.append(len(SeqDictTemp[SeqName]))
		else:
			SSeqName = SSeqName[0].split("_")
			#The Beta names do not need to be shortened:
			if SSeqName[0] == "RefBv":
				NewName = "B_"+SeqName
				NameDictTemp[SeqName] = NewName
				NumBeta += 1
				SeqLenBeta.append(len(SeqDictTemp[SeqName]))				
			else:
				#Give the Portulacineae transcriptomes new transcriptome
				#names that also tell the species.
				TransName = TaxaDict[SSeqName[0]]
				try:
					NewName = TransName+"_"+SSeqName[1]+"_"+SSeqName[2]
				except IndexError:
					NewName = TransName+"_"+SSeqName[1]
				if CatDict[TransName] == 0:
					#plant is Portulacineae, but not Cactaceae
					NewName = "P_"+NewName
					NumPort += 1
					SeqLenPort.append(len(SeqDictTemp[SeqName]))
				elif CatDict[TransName] == 1:
					#plant is Cactaceae
					NewName = "C_"+NewName
					NumCact += 1
					SeqLenCact.append(len(SeqDictTemp[SeqName]))
				elif CatDict[TransName] == 2:
					#plant is outside of the Portulacineae
					NewName = "Y_"+NewName
					NumCary += 1
					SeqLenCary.append(len(SeqDictTemp[SeqName]))
				NameDictTemp[SeqName] = NewName
	#writing the stats for the group to the dictionary
	NumTot = NumMoll + NumBeta + NumPort+NumCact+NumCary
	GroupStats[GroupNum] = [NumTot, NumMoll, SeqLenMoll, NumBeta, SeqLenBeta, NumPort, SeqLenPort, NumCact, SeqLenCact, NumCary, SeqLenCary]
	#Now dealing with the file for the reference genomes:
	InFileName = InFilePre+"Group_"+GroupNum+".fa"
	InFile = open(InFileName, 'rU')
	#resetting everything
	NumRef = 0
	SeqLenRef = [ ]
	PotNum = 1
	#reading the sequences, and renaming potato sequences (the only ones with long names)
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		if Line[0] == ">":
			SeqName = Line[1:]
			if SeqName[0:4] == "PGSC":
				NewName = "R_"+SeqName.split("_")[0]+"_"+str(PotNum)
				NameDictTemp[SeqName] = NewName
				PotNum += 1
			elif SeqName[0:2] == "AT":
				NewName = "A_"+SeqName
				NameDictTemp[SeqName] = NewName
			else:
				NewName = "R_"+SeqName
				NameDictTemp[SeqName] = NewName
		else:
			Seq = Line
			SeqDictTemp[SeqName] = Seq
			NumRef += 1
			SeqLenRef.append(len(Seq))
	#adding the stats for the reference sequences to the stats for that group
	TempList = GroupStats[GroupNum]
	TempList[0] += NumRef
	TempList.append(NumRef)
	TempList.append(SeqLenRef)
	GroupStats[GroupNum] = TempList
	InFile.close()
	NumFiles += 1
	if (NumMoll >= 1) and (NumPort >= 1) and (NumCact >= 1):
		GroupFate[GroupNum] = "accepted! Molluginaceae: %d, Portulacineae %d, Cactaceae %d" % (NumMoll, NumPort, NumCact)
		NumWorked += 1
		ListAccepted.append(GroupNum)
		#writing the sequences to a file
		OutFileName = OutFilePre+GroupNum+".fa"
		OutFile = open(OutFileName, 'w')
		for SeqName in SeqDictTemp.keys():
			Line = ">"+NameDictTemp[SeqName]+"\n"+SeqDictTemp[SeqName]+"\n"
			OutFile.write(Line)
		OutFile.close()
		GroupFileList.append(OutFileName)
		#writing the new names for the sequences to a file
		OutFileName = OutFilePre+GroupNum+"names.txt"
		OutFile = open(OutFileName, 'w')
		for SeqName in NameDictTemp.keys():
			Line = NameDictTemp[SeqName]+"\t"+SeqName+"\n"
			OutFile.write(Line)
		OutFile.close()
	else:
		GroupFate[GroupNum] = "rejected! Molluginaceae: %d, Portulacineae %d, Cactaceae %d" % (NumMoll, NumPort, NumCact)
		NumFailed += 1
	NumFilesOut += 1

if (NumFailed + NumWorked != NumFilesOut):
	print ("Error, %d files read, but only %d of them categorized.\n" % (NumFilesOut, NumFailed+NumWorked))

print("Sequences were read from %d pairs of files; %d of these were accepted and %d were rejected.\n" % (NumFiles, NumWorked, NumFailed))


#writing the stats to another file
#first preparing what needs to be written
StatsWrite = [ ]
Head = "Group\tTotal_Seqs\tNum_Mollugo_Seqs\tLength_Mollugo_Seqs\tNum_Beta_Seqs\tLength_Beta_Seqs\tNum_Portulacineae_non_Cactaceae_Seqs\
\tLength_Portulacineae_non_Cactaceae_Seqs\tNum_Cactaceae_Seqs\tLen_Cactaceae_Seqs\tNum_other_Caryphyllales_seqs\
\tLen_other_Caryophyllales_Seqs\tNum_Reference_Seqs\tLength_Reference_Seqs\tGroup_Fate\n"
StatsWrite.append(Head)
for GroupNum in GroupStats.keys():
	Line = GroupNum
	for Item in GroupStats[GroupNum]:
		Line += "\t"+str(Item)
	Line += "\t"+GroupFate[GroupNum]
	Line += "\n"
	StatsWrite.append(Line)

OutFileName = OutFilePre+"Stats.txt"
OutFile = open(OutFileName, 'w')
for Line in StatsWrite:
	OutFile.write(Line)
OutFile.close()

print("Statistics from %d groups were written to the file %s.\n" % (len(GroupStats.keys()),OutFileName))


#writing the script to build trees with these files:

OutScript = [ ]
RList = [ ]
for GroupNum in ListAccepted:
	#Line = "muscle3.8.31_i86linux64 -in "+OutFilePre+GroupNum+".fa -out "+OutFilePre+GroupNum+"a.fa\n"
	Line = "muscle -in nGroup_"+GroupNum+".fa -out nGroup_"+GroupNum+"a.fa\n"
	OutScript.append(Line)
	Line = OutFilePre+GroupNum+"a.fa\t"+GroupNum+"\n"
	RList.append(Line)
OutFileName = OutFilePre+"RList.txt"
OutFile = open(OutFileName, 'w')
for Line in RList:
	OutFile.write(Line)
OutFile.close()
print ("A list of the alignment files for R is in file %s.\n" % (OutFileName))

OutFileName = OutFilePre+"script"+FileNum+".sh"
OutFile = open(OutFileName, 'w')
for Line in OutScript:
	OutFile.write(Line)
OutFile.close()
print ("A script for aligning these sequences is in file %s.\n" % (OutFileName))

