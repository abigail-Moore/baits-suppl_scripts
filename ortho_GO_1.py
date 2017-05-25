#! /usr/bin/env python

#ortho_GO.py version 1.0 25 Feb. 2014 Abby Moore
#This program is supposed to read the list of groups of orthologous genes made by gene_ortho_2.py and
#find GO terms for those groups, based on the GO terms assigned to the Arabidopsis members.
#It would be nice if it could flag the GO terms to say if they belonged to all AT genes
#in the group, or if they were only specific to one or two.

#The file with GO terms looks like this:
'''
Gene stable ID	Transcript stable ID	GO term accession	GO term name	GO term evidence code	GO domain	Gene name	GO term definition
AT3G18710	AT3G18710.1	GO:0000151	ubiquitin ligase complex	IEA	cellular_component	PUB29	"A protein complex that includes a ubiquitin-protein ligase and other proteins that may confer substrate specificity on the complex." [GOC:jh2, PMID:9529603]
AT3G18710	AT3G18710.1	GO:0003674	molecular_function		molecular_function	PUB29	"Elemental activities, such as catalysis or binding, describing the actions of a gene product at the molecular level. A given gene product may exhibit one or more molecular functions." [GOC:go_curators]
'''
#At this point, we are just interested in columns [0] (AT gene) and [2]
#(GO term, once it is stripped of the GO: part).

#The file with groups (written by gene_ortho_2.py) looks like this:
'''
1	AT5G09800,AT3G18710,AT5G64660	GLYMA06G15630,GLYMA05G32310,GLYMA09G03520,GLYMA08G15580	OS04G0418500,OS05G0439400,OS02G0540700	POPTR_0005s05880,POPTR_0007s03730	PGSC0003DMG400025560	VIT_04s0044g00870
2	AT4G25880,AT3G20250	GLYMA15G17680,GLYMA17G06830,GLYMA13G00670,GLYMA10G28211,GLYMA09G06460	OS08G0519800,OS09G0497100	POPTR_0022s00840,POPTR_0008s00490,POPTR_0006s17870	PGSC0003DMG400009166	VIT_05s0062g00010
'''
#At this point, we are just interested in columns [0] (group number) and
#[1] (AT genese in that group).

import sys #We want to be able to write error messages to the screen.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels.

#The input and output files:
#InFileGroups = "/home/abby/transcriptomes/ortho_groups_head.txt"
InFileGroups = "/home/abby/transcriptomes/orthologous_groups.txt"
#InFileGO = "/home/abby/transcriptomes/GO_terms_head.txt"
InFileGO = "/home/abby/transcriptomes/GO_terms_arabidopsis.plants.ensemble.org.txt"
OutFileGroupsGO = "/home/abby/transcriptomes/groups_to_GO.txt"
OutFileGOGroups = "/home/abby/transcriptomes/GO_to_groups.txt"

#The dictionaries we will need to fill in:
GODict = defaultdict(list)
GroupDict = defaultdict(dict)
GOGroupDict = defaultdict(dict)

#Reading the GO terms to a dictionary where the keys are the AT genes and the values
#are the GO terms.
InFile = open(InFileGO, 'rU')
LineNum = 0
for Line in InFile:
	if LineNum != 0:
		Line = Line.strip('\n').strip('\r').split('\t')
		if Line[2] != "":
			GOTerm = Line[2].split(':')
			GOTerm = GOTerm[1]
			ATGene = Line[0]
			GODict[ATGene].append(GOTerm)
	LineNum += 1
InFile.close()

print ("%d lines were read from the infile %s." % (LineNum, InFileGO))


#Make a default dictionary of groups with the AT gene names as one entry in each group's dictionary:
InFile = open(InFileGroups, 'rU')
LineNum = 0
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	GroupNum = Line[0]
	ATGenes = Line[1].split(',')
	GroupDict[GroupNum]['ATG'] = defaultdict(list)
	GroupDict[GroupNum]['ATG'] = ATGenes
	LineNum += 1
InFile.close()
	
print ("%d lines were read from the infile %s." % (LineNum, InFileGroups))


#Now to join the GO terms with the group numbers:
for GroupNum in GroupDict.keys():
	GroupDict[GroupNum]['universalGO'] = defaultdict(list)
	#If there is only one AT gene in that group, the GO terms for the group are the
	#same as the GO terms for that gene.
	if len(GroupDict[GroupNum]['ATG']) == 1:
		GOListTemp = list(set(GODict[GroupDict[GroupNum]['ATG'][0]]))
		GroupDict[GroupNum]['universalGO'] = GOListTemp
	#If there are multiple genes in the group, it is more complicated.
	elif len(GroupDict[GroupNum]['ATG']) >> 1:
		#first make the list of GO terms that belong to all of the genes in that group:
		GOSetList = [ ]
		for ATGene in GroupDict[GroupNum]['ATG']:
			GOSetList.append(set(GODict[ATGene]))
		GOListTemp = list(set.intersection(*GOSetList))
		GroupDict[GroupNum]['universalGO'] = GOListTemp
		#then make the list of unique/semi-unique GO terms for each gene:
		for ATGene in GroupDict[GroupNum]['ATG']:
			ATGO = GODict[ATGene]
			ATGOListTemp = [ ]
			for GOTerm in ATGO:
				if (GOTerm in GOListTemp) == False:
					ATGOListTemp.append(GOTerm)
			GroupDict[GroupNum][ATGene] = defaultdict(list)
			GroupDict[GroupNum][ATGene] = ATGOListTemp
	else:
		print ("There is an error with group %s!!" % (GroupNum))

		
#Now to make the output file linking the groups to their GO terms:

FileList = [ ]
Line = "Group number\tArabidopsis genes in group\tUniversal GO terms\tGO terms specific to each gene\n"
FileList.append(Line)
LineNum = 0
for GroupNum in GroupDict.keys():
	Line = GroupNum+"\t"
	for ATGene in GroupDict[GroupNum]['ATG']:
		Line += ATGene+","
		#getting rid of the last comma:
	Line = Line[:-1]
	Line += '\t'
	for GOTerm in GroupDict[GroupNum]['universalGO']:
		Line += GOTerm+","
	Line = Line[:-1]
	Line += '\t'
	if len(GroupDict[GroupNum]['ATG']) >> 1:
		for ATGene in GroupDict[GroupNum]['ATG']:
			Line += ATGene+":"
			if GroupDict[GroupNum][ATGene] == [ ]:
				Line += "no unique GO terms\t"
			else:
				for GOTerm in GroupDict[GroupNum][ATGene]:
					Line += GOTerm+","
				Line = Line[:-1]
				Line += "\t"
		Line = Line[:-1]
	Line += "\n"
	FileList.append(Line)
	LineNum += 1

OutFile = open(OutFileGroupsGO, 'w')
for Line in FileList:
	OutFile.write(Line)
OutFile.close()

print ("Information for %d groups was written to the outfile %s." % (LineNum, OutFileGroupsGO))

#Now to join the group numbers with the GO terms (same thing, just reversed):

for GroupNum in GroupDict.keys():
	for GOTerm in GroupDict[GroupNum]['universalGO']:
		try:
			#All of this difficulty is because I cannot figure out how to append
			#directly to a defaultdict(list).
			ListTemp = GOGroupDict[GOTerm]['universal']
			ListTemp.append(GroupNum)
			GOGroupDict[GOTerm]['universal'] = ListTemp
		except KeyError:
			GOGroupDict[GOTerm]['universal'] = defaultdict(list)
			GOGroupDict[GOTerm]['universal'] = [GroupNum]
	if len(GroupDict[GroupNum]['ATG']) >> 1:
		for ATGene in GroupDict[GroupNum]['ATG']:
			for GOTerm in GroupDict[GroupNum][ATGene]:
				try:
					ListTemp = GOGroupDict[GOTerm]['specific']
					ListTemp.append(GroupNum)
					GOGroupDict[GOTerm]['specific'] = ListTemp
				except KeyError:
					GOGroupDict[GOTerm]['specific'] = defaultdict(list)
					GOGroupDict[GOTerm]['specific'] = [GroupNum]
#Now the duplicate group numbers need to be removed from the list.
for GOTerm in GOGroupDict.keys():
	#The specific terms could well be present multiple times, so duplicates need to be removed
	#before sorting.
	try:
		ListTemp = GOGroupDict[GOTerm]['specific']
		ListTemp = list(set(ListTemp))
		GOGroupDict[GOTerm]['specific'] = ListTemp
	except KeyError:
		'something'

FileList = [ ]
Line = "GO Term\tAll Have\tSome Have\n"
FileList.append(Line)
LineNum = 0
KeysSorted = GOGroupDict.keys()
KeysSorted.sort()
for GOTerm in KeysSorted:
	Line = GOTerm+"\t"
	try: 
		if len(GOGroupDict[GOTerm]['universal']) != 0:
			for GroupNum in GOGroupDict[GOTerm]['universal']:
				Line += GroupNum+","
			Line = Line[:-1]
			Line += "\t"
	except KeyError:
		Line += "none\t"
	try: 
		if len(GOGroupDict[GOTerm]['specific']) != 0:
			for GroupNum in GOGroupDict[GOTerm]['specific']:
				Line += GroupNum+","
			Line = Line[:-1]
	except KeyError:
		Line += "none\t"
	Line += "\n"
	FileList.append(Line)
	LineNum += 1
	
OutFile = open(OutFileGOGroups, 'w')
for Line in FileList:
	OutFile.write(Line)
OutFile.close()

print ("Information for %d GO terms was written to the outfile %s." % (LineNum, OutFileGOGroups))
