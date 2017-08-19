#! /usr/bin/env python

#talign_groups.py version 1.0 1 May 2014 Abby Moore
#This script is supposed to look at the final alignments to determine whether or not they
#have Arabidopsis sequences (and thus whether or not GO terms can be determined).  Then,
#for the alignments that have Arabidopsis sequences, it assigns them GO terms in three
#categories of certainty:
#1: all of the genes in the group of homologous sequences (of which the aligned sequences are
#in most cases a subset) have that GO term
#2: all of the Arabidopsis sequences in this alignment (but not all of them in the group
#as a whole) have that GO term
#3: one or more, but not all, of the Arabidopsis sequences in this alignment have that
#GO term

#File Formats:
#GO_to_groups.txt (tab delimited, all on one line, first line is header):
#although I am not sure we actually need this file??
#0000001: GO term [0]
#13861,8070,13067,4159,22240,9350,24020: groups in which all sequences have GO term, separated by commas [1]
#15752,290: groups in which some of the sequences have GO term, separated by commas [2]

#groups_to_GO.txt (tab delimited, all on one line, first line is header):
#11542: group number [0]
#AT2G40000,AT3G55840: Arabidopsis genes in group, separated by commas [1]
#0006950,0000084....: GO terms all Arabidopsis genes in group have, separated by commas [2]
#AT2G40000:0006979,0009646....: GO terms this locus has that all of them don't have, separated by commas [3]
#AT3G55840:0002679,0007165....: same for second locus [4], etc.

from collections import defaultdict #We want to be able to create dictionaries with multiple levels
import re #We will need regular expressions to read the group names.

#Input files:
#ALListFile = "/home/abby/transcriptomes/finalgroupshead.txt"
#ALListFile = "/home/abby/transcriptomes/finalgroupshead20.txt"
ALListFile = "/home/abby/transcriptomes/finalgroups.txt"
ALFilePre = "/home/abby/transcriptomes/finalgroups/"
GOGroupsFile = "/home/abby/transcriptomes/GO_to_groups.txt"
GroupsGOFile = "/home/abby/transcriptomes/groups_to_GO.txt"
#GroupsGOFile = "/home/abby/transcriptomes/headGtoGO.txt"

#Output files:
ALFileOut = "/home/abby/transcriptomes/finalgroups/ALkey.txt"
GOFileOut = "/home/abby/transcriptomes/finalgroups/GOkey.txt"
ALSeqsList = "/home/abby/transcriptomes/finalgroups/SeqsperAL.txt"
GOSpreadList = "/home/abby/transcriptomes/finalgroups/GOspread.txt"

#Lists and dictionaries we will fill out:
ALList = [ ] #List of all alignment files
GroupGODict = defaultdict(dict) #Dictionary written from the groups_to_GO file
ALDict = defaultdict(dict) #Dictionary of all of the information sorted by alignment file
GOALDict = defaultdict(dict) #Dictionary of all of the information sorted by GO term.

GroupRe = r"(\d+)_\w+"
GroupSub = r"\1"
ALRe = r"(\d+_)\w(_\w+\.fa)"
ALSub = r"\1B\2"

#Reading the list of files and saving them as a list:
InFile = open(ALListFile, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	ALList.append(Line)
InFile.close()
print ("Infile %s read." % (ALListFile))

#Looking through the alignment files to see if they have Arabidopsis genes in them,
#and making a list of Portullugo genes:
MOListTaken = [ ]
for FileName in ALList:
	InFileName = ALFilePre + FileName
	ALName = FileName[7:-4]
	ATListTemp = [ ]
	MOListTemp = [ ]
	CactListTemp = [ ]
	PortListTemp = [ ]
	AlignDictTemp = { }
	InFile = open(InFileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\n').strip('\r')
		if Line[0] == ">":
			SeqName = Line
		else:
			try:
				AlignDictTemp[SeqName] += Line
			except KeyError:
				AlignDictTemp[SeqName] = Line
	InFile.close()
	for SeqName in AlignDictTemp.keys():
		SeqLen = 0
		for Char in AlignDictTemp[SeqName]:
			if Char != "-":
				SeqLen += 1
		if SeqName[0:2] == ">A":
			ATListTemp.append(SeqName[3:-2])
		elif SeqName[0:2] == ">M":
			Entry = SeqName[3:]+":"+str(SeqLen)
			MOListTemp.append(Entry)
		elif SeqName[0:2] == ">C":
			Entry = SeqName[3:]+":"+str(SeqLen)
			CactListTemp.append(Entry)
		elif SeqName[0:2] == ">P":
			Entry = SeqName[3:]+":"+str(SeqLen)
			PortListTemp.append(Entry)
	if ((MOListTemp[0] in MOListTaken) == False):
		for MOTerm in MOListTemp:
			MOListTaken.append(MOTerm)
		if ATListTemp == [ ]:
			ALDict[ALName]['AT'] = "none present"
		else:
			ATListTemp = list(set(ATListTemp))
			ALDict[ALName]['AT'] = ATListTemp
		ALDict[ALName]['MO'] = MOListTemp
		ALDict[ALName]['Cact'] = CactListTemp
		ALDict[ALName]['Port'] = PortListTemp
	#else:
	#	print ("%s rejected!" % (ALName))

#Reading the groups_to_GO file and making a dictionary
InFile = open(GroupsGOFile, 'rU')
LineNum = 1
for Line in InFile:
	if LineNum > 1:
		Line = Line.strip('\n').strip('\r').split('\t')
		GroupNum = Line[0]
		NumAT = len(Line[1].split(','))
		if NumAT == 1:
			GroupGODict[GroupNum]['Univ'] = Line[2].split(',')
		if NumAT > 1:
			#if there are some universal GO terms
			if len(Line) == 3+NumAT:
				#add the universal GO terms to the dictionary
				GroupGODict[GroupNum]['Univ'] = Line[2].split(',')
				#look through the unique GO terms and give each AT gene
				#its own entry
				for Num in range(3,3+NumAT):
					SubList = Line[Num].split(':')
					ATGene = SubList[0]
					ATGOList = SubList[1].split(',')
					GroupGODict[GroupNum][ATGene] = ATGOList
			#if there are no universal GO terms (and the line is therefore
			#shorter)
			else:
				GroupGODict[GroupNum]['Univ'] = [ ]
				for Num in range(2,2+NumAT):
					SubList = Line[Num].split(':')
					ATGene = SubList[0]
					ATGOList = SubList[1].split(',')
					GroupGODict[GroupNum][ATGene] = ATGOList
	LineNum += 1
InFile.close()
print("Infile %s read." % (GroupsGOFile))

#matching the GO terms to the alignments
for ALName in ALDict:
	#only necessary of this alignment actually has Arabidopsis genes
	if ALDict[ALName]['AT'] != "none present":
		#finding out which group the alignment belongs to
		GroupNum = re.sub(GroupRe, GroupSub, ALName)
		GOListTemp = [ ]
		GOListAll = [ ]
		#adding the GO terms that belong to all members of that group
		if GroupGODict[GroupNum]['Univ'] == [ ]:
			ALDict[ALName]['Univ'] = ["none present"]
		else:
			ALDict[ALName]['Univ'] = GroupGODict[GroupNum]['Univ']
		#finding the GO terms that belong to some members of the group
		for ATGene in ALDict[ALName]['AT']:
			try:
				if GroupGODict[GroupNum][ATGene] != ["no unique GO terms"]:
					GOListTemp += GroupGODict[GroupNum][ATGene]
			except KeyError:
				"No GO Terms"
		GOListSet = list(set(GOListTemp))
		#determining of any of these belong to all Arabidposis sequences that are
		#part of that alignment
		for GOTerm in GOListSet:
			if GOListTemp.count(GOTerm) == len(ALDict[ALName]['AT']):
				GOListAll.append(GOTerm)
				GOListSet.remove(GOTerm)
		ALDict[ALName]['All'] = GOListAll
		ALDict[ALName]['Some'] = GOListSet

#Now reversing that dictionary to make the dictionary we really want
for ALName in ALDict:
	if ALDict[ALName]['AT'] != "none present":
		try:
			for GOTerm in ALDict[ALName]['Univ']:
				try:
					GOALDict[GOTerm]['Univ'].append(ALName)
				except KeyError:
					GOALDict[GOTerm]['Univ'] = [ALName]
		except KeyError:
			"no universal GO terms"
		for GOTerm in ALDict[ALName]['Some']:
			try:
				GOALDict[GOTerm]['Some'].append(ALName)
			except KeyError:
				GOALDict[GOTerm]['Some'] = [ALName]
		for GOTerm in ALDict[ALName]['All']:
			try:
				GOALDict[GOTerm]['All'].append(ALName)
			except KeyError:
				GOALDict[GOTerm]['All'] = [ALName]

#If there are lists of GO terms that are especially desired or list of GO term functions,
#they could be added to the GOALDict now.

#printing the dictionary by alignment name
ALFileList = [ ]
ALNameList = ALDict.keys()
ALNameList.sort()
for ALName in ALNameList:
	Line = ALName + "\t"
	for MOEntry in ALDict[ALName]['MO']:
		Line += MOEntry+","
	Line = Line[:-1]+"\t"
	for CactEntry in ALDict[ALName]['Cact']:
		Line += CactEntry+","
	Line = Line[:-1]+"\t"
	for PortEntry in ALDict[ALName]['Port']:
		Line += PortEntry+","
	Line = Line[:-1]+"\t"
	if ALDict[ALName]['AT'] == "none present":
		Line += "none present\n"
	elif ALDict[ALName]['AT'] != "none present":
		for ATName in ALDict[ALName]['AT']:
			Line += ATName+","
		Line = Line[:-1]+"\t"
		for GOTerm in ALDict[ALName]['Univ']:
			Line += GOTerm+","
		Line = Line[:-1]+"\t"
		if ALDict[ALName]['All'] == "":
			Line += "\t"
		else:
			for GOTerm in ALDict[ALName]['All']:
				Line += GOTerm+","
				Line = Line[:-1]+"\t"
		if ALDict[ALName]['Some'] == "":
			Line += "\t"
		else:
			for GOTerm in ALDict[ALName]['Some']:
				Line += GOTerm+","
				Line = Line[:-1]+"\t"
		Line += "\n"
	ALFileList.append(Line)

OutFile = open(ALFileOut, 'w')
for Line in ALFileList:
	OutFile.write(Line)
OutFile.close()
print("Outfile %s written." % (ALFileOut))

#printing the dictionary by GO term
GOFileList = [ ]
GOTermList = GOALDict.keys()
GOTermList.sort()
GOTermList.remove("")
for GOTerm in GOTermList:
	Line = GOTerm + "\t"
	try:
		for ALName in GOALDict[GOTerm]['Univ']:
			Line += ALName + ","
		Line = Line[:-1]+"\t"
	except KeyError:
		Line += "\t"
	try:
		for ALName in GOALDict[GOTerm]['All']:
			Line += ALName+","
		Line = Line[:-1]+"\t"
	except KeyError:
		Line += "\t"
	try:
		for ALName in GOALDict[GOTerm]['Some']:
			Line += ALName+","
		Line = Line[:-1]+"\t"
	except KeyError:
		Line += "\t"
	Line += "\n"
	GOFileList.append(Line)

OutFile = open(GOFileOut, 'w')
for Line in GOFileList:
	OutFile.write(Line)
OutFile.close()
print("Outfile %s written." % (GOFileOut))

#printing the file that just lists the sequences by alignment (with their lengths),
#so this can be modified when choosing sequences:
ALSeqsFile = [ ]
ALNameList = ALDict.keys()
ALNameList.sort()
for ALName in ALNameList:
	for MOEntry in ALDict[ALName]['MO']:
		Line = ALName+"\t"+MOEntry+"\n"
		ALSeqsFile.append(Line)
	for CactEntry in ALDict[ALName]['Cact']:
		Line = ALName+"\t"+CactEntry+"\n"
		ALSeqsFile.append(Line)
	for PortEntry in ALDict[ALName]['Port']:
		Line = ALName+"\t"+PortEntry+"\n"
		ALSeqsFile.append(Line)
OutFile = open(ALSeqsList, 'w')
for Line in ALSeqsFile:
	OutFile.write(Line)
OutFile.close()
print("Outfile %s written." % (ALSeqsList))

#printing a spread out GO term file:
GOSpreadFile = [ ]
for GOTerm in GOTermList:
	try:
		for ALName in GOALDict[GOTerm]['Univ']:
			Line = GOTerm+"\t"+ALName+"\tUniv\t"
			for otherTerm in ALDict[ALName]['Univ']:
				Line += otherTerm+","
			Line = Line[:-1]+"\n"
			GOSpreadFile.append(Line)
	except KeyError:
		"do nothing"
	try:
		for ALName in GOALDict[GOTerm]['All']:
			Line = GOTerm+"\t"+ALName+"\tAll\t"
			for otherTerm in ALDict[ALName]['Univ']:
				Line += otherTerm+","
			Line = Line[:-1]+"\n"
			GOSpreadFile.append(Line)
	except KeyError:
		"do nothing"
	try:
		for ALName in GOALDict[GOTerm]['Some']:
			Line = GOTerm+"\t"+ALName+"\tSome\t"
			for otherTerm in ALDict[ALName]['Univ']:
				Line += otherTerm+","
			Line = Line[:-1]+"\n"
			GOSpreadFile.append(Line)
	except KeyError:
		"do nothing"


OutFile = open(GOSpreadList, 'w')
for Line in GOSpreadFile:
	OutFile.write(Line)
OutFile.close()
print("Outfile %s written.\n" % (GOSpreadList))
