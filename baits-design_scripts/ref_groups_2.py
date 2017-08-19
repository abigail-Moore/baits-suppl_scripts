#! /usr/bin/env python

#ref_groups.py version 2.0 6 March 2014 Abby Moore
#This file reads the transcriptome files from Ensembl and the list of sequences in each
#group of orthologs and finds the sequences that belong to the different groups.
#The sequences of each group are saved in separate files and also in a combined file, with the
#sequences in that combined file labeled as to which group they are in.

import re #We will need the regular expressions module
import sys #We want to be able to output error messages to the screen and get input from the command line.
import gzip #We want to be able to open zipped files.
from collections import defaultdict #We want to make a dictionary with multiple levels

#The form of the infile for the groups (and thus the required formats for the names):
'''
1	AT5G09800,AT3G18710,AT5G64660	GLYMA06G15630,GLYMA05G32310,GLYMA09G03520,GLYMA08G15580	OS04G0418500,OS05G0439400,OS02G0540700	POPTR_0005s05880,POPTR_0007s03730	PGSC0003DMG400025560	VIT_04s0044g00870
2	AT4G25880,AT3G20250	GLYMA15G17680,GLYMA17G06830,GLYMA13G00670,GLYMA10G28211,GLYMA09G06460	OS08G0519800,OS09G0497100	POPTR_0022s00840,POPTR_0008s00490,POPTR_0006s17870	PGSC0003DMG400009166	VIT_05s0062g00010
3	AT1G71695	GLYMA16G27900,GLYMA10G36690,GLYMA16G27880,GLYMA02G08780,GLYMA16G27890				VIT_18s0072g00160
'''

#NOTE that when I change to the full files, I will have to change "open" to "gzip.open"!!!

#The infiles:

#The infiles for the transcriptomes:
InFileAT = "/home/abby/transcriptomes/reference_transcriptomes/Arabidopsis_thaliana.TAIR10.21.cdna.all.fa.gz"
#InFileAT = "/home/abby/transcriptomes/reference_transcriptomes/Arabidopsis_thaliana_completecds.fa.gz"
#InFileAT = "/home/abby/transcriptomes/reference_transcriptomes/Arabidopsis_thaliana_head.fa"
InFileGM = "/home/abby/transcriptomes/reference_transcriptomes/Glycine_max.V1.0.21.cdna.all.fa.gz"
#InFileGM = "/home/abby/transcriptomes/reference_transcriptomes/Glycine_max.fa.gz"
#InFileGM = "/home/abby/transcriptomes/reference_transcriptomes/Glycine_max_head.fa"
InFileOS = "/home/abby/transcriptomes/reference_transcriptomes/Oryza_sativa.IRGSP-1.0.21.cdna.all.fa.gz"
#InFileOS = "/home/abby/transcriptomes/reference_transcriptomes/Osativa_head.fa"
InFilePT = "/home/abby/transcriptomes/reference_transcriptomes/Populus_trichocarpa.JGI2.0.21.cdna.all.fa.gz"
#InFilePT = "/home/abby/transcriptomes/reference_transcriptomes/Ptrichocarpa_210_cds.fa.gz"
#InFilePT = "/home/abby/transcriptomes/reference_transcriptomes/Ptrichocarpa_210_head.fa"
#InFilePTsyn = "/home/abby/transcriptomes/reference_transcriptomes/Populus_trichocarpa_synonyms.txt"
InFileST = "/home/abby/transcriptomes/reference_transcriptomes/Solanum_tuberosum.3.0.21.cdna.all.fa.gz"
#InFileST = "/home/abby/transcriptomes/reference_transcriptomes/Stuberosum_head.fa"
InFileVV = "/home/abby/transcriptomes/reference_transcriptomes/Vitis_vinifera.IGGP_12x.21.cdna.all.fa.gz"
#InFileVV = "/home/abby/transcriptomes/reference_transcriptomes/Vitis_vinifera_head.fa"

#The infile for the groups:
InFileGroups = "/home/abby/transcriptomes/orthologous_groups.txt"
#InFileGroups = "/home/abby/transcriptomes/ortho_groups_head.txt"

#The outfile(s?):
OutFilePre = "/home/abby/transcriptomes/groups/Group_"
OutFilePost = ".fa"

#The dictionaries to be made with the infiles:
ATDict = defaultdict(dict)
GMDict = defaultdict(dict)
OSDict = defaultdict(dict)
PTDict = defaultdict(dict)
STDict = defaultdict(dict)
VVDict = defaultdict(dict)
GroupDict = defaultdict(dict)

#Reading the Arabidopsis thaliana sequences.  The names of the sequences are in the same format
#as the names in the list of groups.
SeqNum = 0
Seq = ""
InFile = gzip.open(InFileAT)
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	if Line[0] == ">":
		#Before starting with the next sequence, add the previous one to the dictionary
		if Seq != "":
			ATDict[GeneName][SeqName] = Seq
			SeqNum += 1
		#Now we can save the name of the next sequence:
		SeqName = Line.strip('>').split(' ')[0]
		GeneName = SeqName.split('.')[0]
		Seq = ""
	else:
		#The sequences are spread over more than one line
		Seq += Line
ATDict[GeneName][SeqName] = Seq
SeqNum += 1
InFile.close()
print ("%d sequences were read from the infile %s.\n" % (SeqNum, InFileAT))

#Reading the Glycine max sequences.
SeqNum = 0
Seq = ""
InFile = gzip.open(InFileGM)
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	if Line[0] == ">":
		if Seq != "":
			GMDict[GeneName][SeqName] = Seq
			SeqNum += 1
		SeqName = Line.strip('>').split(' ')[0]
		GeneName = SeqName.split('.')[0]
		Seq = ""
	else:
		Seq += Line
GMDict[GeneName][SeqName] = Seq
SeqNum += 1
InFile.close()
print ("%d sequences were read from the infile %s.\n" % (SeqNum, InFileGM))

#Reading the Oryza sativa sequences.
#format:
#>OS01T0100100-01 cdna:known chromosome:IRGSP-1.0:1:2983:10815:1 gene:OS01G0100100 transcript:OS01T0100100-01
SeqNum = 0
Seq = ""
InFile = gzip.open(InFileOS)
#InFile = open(InFileOS, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	if Line[0] == ">":
		if Seq != "":
			OSDict[GeneName][SeqName] = Seq
			SeqNum += 1
		Line = Line.strip('>').split(' ')
		GeneName = Line[3].split(':')[1]
		SeqName = GeneName+"_"+Line[4].split(':')[1]
		Seq = ""
	else:
		Seq += Line
OSDict[GeneName][SeqName] = Seq
SeqNum += 1
InFile.close()
print ("%d sequences were read from the infile %s.\n" % (SeqNum, InFileOS))

#Reading the Populus trichocarpa sequences.
SeqNum = 0
Seq = ""
InFile = gzip.open(InFilePT)
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	if Line[0] == ">":
		if Seq != "":
			PTDict[GeneName][SeqName] = Seq
			SeqNum += 1
		SeqName = Line.strip('>').split(' ')[0]
		GeneName = SeqName.split('.')[0]
		Seq = ""
	else:
		Seq += Line
PTDict[GeneName][SeqName] = Seq
SeqNum += 1
InFile.close()
print ("%d sequences were read from the infile %s.\n" % (SeqNum, InFilePT))

#Reading the Solanum tuberosum sequences.
#format:
#>PGSC0003DMT400039136 cdna:novel chromosome:3.0:1:152322:153489:-1 gene:PGSC0003DMG400015133 transcript:PGSC0003DMT400039136 description:"Defensin"
SeqNum = 0
Seq = ""
InFile = gzip.open(InFileST)
#InFile = open(InFileST, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	if Line[0] == ">":
		if Seq != "":
			STDict[GeneName][SeqName] = Seq
			SeqNum += 1
		Line = Line.strip('>').split(' ')
		GeneName = Line[3].split(':')[1]
		SeqName = GeneName+"_"+Line[4].split(':')[1]
		Seq = ""
	else:
		Seq += Line
STDict[GeneName][SeqName] = Seq
SeqNum += 1
InFile.close()
print ("%d sequences were read from the infile %s.\n" % (SeqNum, InFileST))

#Reading the Vitis vinifera sequences.  The names of the sequences are in the same format
#as the names in the list of groups.
SeqNum = 0
Seq = ""
#InFile = open(InFileVV, 'rU')
InFile = gzip.open(InFileVV)
for Line in InFile:
	Line = Line.strip('\r').strip('\n')
	if Line[0] == ">":
		if Seq != "":
			VVDict[GeneName][SeqName] = Seq
			SeqNum += 1
		SeqName = Line.strip('>').split(' ')[0]
		GeneName = SeqName.split('.')[0]
		Seq = ""
	else:
		Seq += Line
VVDict[GeneName][SeqName] = Seq
SeqNum += 1
InFile.close()
print ("%d sequences were read from the infile %s.\n" % (SeqNum, InFileVV))

#Reading the file with information on the various groups and making a dictionary with the sequences:
NumGroups = 0
NumMissing = 0
NumGroupsMissing = 0
GroupList = [ ]
InFile = open(InFileGroups, 'rU')
for Line in InFile:
	Line = Line.strip('\r').strip('\n').split('\t')
	GroupNum = Line[0]
	GroupList.append(GroupNum)
	if Line[1] != "":
		GeneList = Line[1].split(',')
		for GeneName in GeneList:
			if ATDict[GeneName].keys() == [ ]:
				print ("Error, gene %s has no sequence!!\n" % (GeneName))
				NumMissing += 1
			for SeqName in ATDict[GeneName].keys():
				GroupDict[GroupNum][SeqName] = ATDict[GeneName][SeqName]
	if Line[2] != "":
		GeneList = Line[2].split(',')
		for GeneName in GeneList:
			if GMDict[GeneName].keys() == [ ]:
				print ("Error, gene %s has no sequence!!\n" % (GeneName))
				NumMissing += 1
			for SeqName in GMDict[GeneName].keys():
				GroupDict[GroupNum][SeqName] = GMDict[GeneName][SeqName]
	if Line[3] != "":
		GeneList = Line[3].split(',')
		for GeneName in GeneList:
			if OSDict[GeneName].keys() == [ ]:
				print ("Error, gene %s has no sequence!!\n" % (GeneName))
				NumMissing += 1
			for SeqName in OSDict[GeneName].keys():
				GroupDict[GroupNum][SeqName] = OSDict[GeneName][SeqName]
	if Line[4] != "":
		GeneList = Line[4].split(',')
		for GeneName in GeneList:
			if PTDict[GeneName].keys() == [ ]:
				print ("Error, gene %s has no sequence!!\n" % (GeneName))
				NumMissing += 1
			for SeqName in PTDict[GeneName].keys():
				GroupDict[GroupNum][SeqName] = PTDict[GeneName][SeqName]
	if Line[5] != "":
		GeneList = Line[5].split(',')
		for GeneName in GeneList:
			if STDict[GeneName].keys() == [ ]:
				print ("Error, gene %s has no sequence!!\n" % (GeneName))
				NumMissing += 1
			for SeqName in STDict[GeneName].keys():
				GroupDict[GroupNum][SeqName] = STDict[GeneName][SeqName]
	if Line[6] != "":
		GeneList = Line[6].split(',')
		for GeneName in GeneList:
			if VVDict[GeneName].keys() == [ ]:
				print ("Error, gene %s has no sequence!!\n" % (GeneName))
				NumMissing += 1
			for SeqName in VVDict[GeneName].keys():
				GroupDict[GroupNum][SeqName] = VVDict[GeneName][SeqName]
	if ((Line[1] == "") and (Line[2] == "") and (Line[3] == "") and (Line[4] == "") and (Line[5] == "") and (Line[6] == "")):
		NumGroupsMissing += 1
	NumGroups += 1
InFile.close()

print ("Information about %d groups of orthologous sequences was read from the file %s.\n" % (NumGroups, InFileGroups))
print ("The dictionary of these sequences is of length %d, after empty groups were excluded.\n" % (len(GroupDict.keys())))

print ("%d genes could not be found in the transcriptomes.\n" % (NumMissing))

#Printing the grouped sequences to a file.  Here I am printing all of the sequences to a single file,
#but am naming the sequences according their groups.

SeqFile = [ ]
NumSeqs = 0
NumGroups = 0
for GroupNum in GroupDict.keys():
	for SeqName in GroupDict[GroupNum].keys():
		Seq = ">"+GroupNum+"_"+SeqName+"\n"+GroupDict[GroupNum][SeqName]+"\n"
		SeqFile.append(Seq)
		NumSeqs += 1
	NumGroups += 1

OutFileName = OutFilePre+"All"+OutFilePost
OutFile = open(OutFileName, "w")
for Seq in SeqFile:
	OutFile.write(Seq)
OutFile.close()

print ("Information for %d sequences in %d groups was written to the file %s.\n" % (NumSeqs, NumGroups, OutFileName))

#Printing the sequences to individual files for each group:

NumSeqs = 0
NumGroups = 0
for GroupNum in GroupDict.keys():
	OutFileName = OutFilePre+GroupNum+OutFilePost
	OutFile = open(OutFileName, 'w')
	for SeqName in GroupDict[GroupNum].keys():
		Seq = ">"+SeqName+"\n"+GroupDict[GroupNum][SeqName]+"\n"
		OutFile.write(Seq)
		NumSeqs += 1
	OutFile.close()
	NumGroups += 1

print ("Information for the same %d sequences in the same %d groups was also " % (NumSeqs, NumGroups))
print ("written to %d separate files, one for each group, \nwith names such as %s.\n" % (NumGroups, OutFileName))	
