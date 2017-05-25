#! /usr/bin/env python

#This program is supposed to read the gene_orthology_modelspecies.plants.ensemble.org.txt file,
#make groups of orthologous genes, and then make separate files for each of the sets of genes.
#gene_ortho_2.py Version 2.0, 24 Feb. 2014, Abby Moore

#The input file looks like this:
'''Gene stable ID	Transcript stable ID	Glycine max Ensembl Gene ID	Homology Type	% Identity	Glycine max % Identity	Oryza sativa Japonica (Rice) Ensembl Gene ID	Homology Type	% Identity	Oryza sativa Japonica (Rice) % Identity	Populus trichocarpa Ensembl Gene ID	Homology Type	% Identity	Populus trichocarpa % Identity	Solanum tuberosum Ensembl Gene ID	Homology Type	% Identity	Solanum tuberosum % Identity	Vitis vinifera Ensembl Gene ID	Homology Type	% Identity	Vitis vinifera % Identity
AT3G18710	AT3G18710.1	GLYMA05G32310	ortholog_many2many	50	49	OS04G0418500	ortholog_many2many	30	27	POPTR_0005s05880	ortholog_many2many	50	51	PGSC0003DMG400025560	ortholog_one2many	42	41	VIT_04s0044g00870	ortholog_one2many	50	50
AT3G18710	AT3G18710.1	GLYMA08G15580	ortholog_many2many	49	49	OS04G0418500	ortholog_many2many	30	27	POPTR_0005s05880	ortholog_many2many	50	51	PGSC0003DMG400025560	ortholog_one2many	42	41	VIT_04s0044g00870	ortholog_one2many	50	50
'''

import sys #We will want to be able to output error messages to the screen.
from collections import defaultdict #We want to be able to make dictionaries with multiple levels.

#The input and output files:
#InFileGeneOrth = "/home/abby/transcriptomes/gene_orthology_head.txt"
InFileGeneOrth = "/home/abby/transcriptomes/gene_orthology_modelspecies.plants.ensemble.org.txt"
OutFileGeneOrth = "/home/abby/transcriptomes/orthologous_groups.txt"

#The lists and dictionaries of we will make from the input file:
GrpDict = defaultdict(dict) #This is a default dict, because each group will have separate entries
#for AT, GM, OS, PT, ST, and VV genes.
GMAT = { } #These are dictionaries that say which AT genes are orthologous with the genes from the other
#model plants.
OSAT = { }
PTAT = { }
STAT = { }
VVAT = { }
ATGM = { }
ATOS = { }
ATPT = { }
ATST = { }
ATVV = { }
#These are the lists of genes.  The non-AT genes have two lists, because the "temp" list
#will have most all genes listed multiple times, and the other list is the corrected version.
ATList = [ ]
GMListtemp = [ ]
GMList = [ ]
OSListtemp = [ ]
OSList = [ ]
PTListtemp = [ ]
PTList = [ ]
STListtemp = [ ]
STList = [ ]
VVListtemp = [ ]
VVList = [ ]
#These are the dictionaries that will have each gene together with its group number.
ATGrpDict = { }
GMGrpDict = { }
OSGrpDict = { }
PTGrpDict = { }
STGrpDict = { }
VVGrpDict = { }

#The lists and dictionaries we will output:
GroupDict = { } #This is the final dictionary of groups of orthologous genes, to be written to
#OutFileGeneOrth.

#The first thing to do is to read the input file and fill out the dictionary and lists.
#What I have done here is a bit clumsy, but, as far as I know, there is not a way to loop over 
#dictionary names.
InFile = open(InFileGeneOrth, 'rU')
LineNum = 0
ATGene0 = ''
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	#This gets rid of the line ending characters and changes the line from a string 
	#(whose elements were originally separated by tabs) into a list.
	if LineNum != 0: #We want to skip the header row
		ATGene = Line[0]
		if ATGene != ATGene0:
			ATList.append(ATGene)
			ATGene0 = ATGene
		GMGene = Line[2]
		if GMGene != "": #If we have a GM gene...
			#add the gene to the list of GM genes
			GMListtemp.append(GMGene)
			#First, deal with the GM-AT dictionary:
			try: #Try to add the AT gene to the existing entry for that GM gene.
				GMAT[GMGene] += ","+ATGene
			except KeyError: #Make a new entry for that GM gene if it doesn't have an entry yet.
				GMAT[GMGene] = ATGene
			#Then deal with the AT-GM dictionary:
			try: #Try to add the GM gene to the existing entry for that AT gene.
				ATGM[ATGene] += ","+GMGene
			except KeyError: #Make a new entry for that AT gene if it doesn't have an entry yet.
				ATGM[ATGene] = GMGene
		OSGene = Line[6]
		if OSGene != "": #And so on for OS, PT, ST, and VV.....
			OSListtemp.append(OSGene)
			try:
				OSAT[OSGene] += ","+ATGene
			except KeyError:
				OSAT[OSGene] = ATGene
			try:
				ATOS[ATGene] += ","+OSGene
			except KeyError:
				ATOS[ATGene] = OSGene
		PTGene = Line[10]
		if PTGene != "":
			PTListtemp.append(PTGene)
			try:
				PTAT[PTGene] += ","+ATGene
			except KeyError:
				PTAT[PTGene] = ATGene
			try:
				ATPT[ATGene] += ","+PTGene
			except KeyError:
				ATPT[ATGene] = PTGene
		STGene = Line[14]
		if STGene != "":
			STListtemp.append(STGene)
			try:
				STAT[STGene] += ","+ATGene
			except KeyError:
				STAT[STGene] = ATGene
			try:
				ATST[ATGene] += ","+STGene
			except KeyError:
				ATST[ATGene] = STGene
		VVGene = Line[18]
		if VVGene != "":
			VVListtemp.append(VVGene)
			try:
				VVAT[VVGene] += ","+ATGene
			except KeyError:
				VVAT[VVGene] = ATGene
			try:
				ATVV[ATGene] += ","+VVGene
			except KeyError:
				ATVV[ATGene] = VVGene		
		#print ("%s AT, %s GM, %s OS, %s PT, %s ST, %s VV.\n" % (ATGene, GMGene, OSGene, PTGene, STGene, VVGene))
	LineNum += 1
InFile.close()

print ("%d lines were read from the file %s.\n" % (LineNum, InFileGeneOrth))
#Making versions of the various (non-AT) lists where each gene is present only once.
#This is almost certainly faster than checking whether a gene is present each time we want to add
#it to the list.

GMList = list(set(GMListtemp))
OSList = list(set(OSListtemp))
PTList = list(set(PTListtemp))
STList = list(set(STListtemp))
VVList = list(set(VVListtemp))

#Now to give everything a group:
GrpNum = 1
for ATGene in ATList: #Going through the list of AT genes.
	try: #First, see if that AT gene already has a group number (because one of its orthologs
		#in another species has already come up).
		ATNum = ATGrpDict[ATGene]
	except KeyError:
		#If not, we start a new group number.
		ATNum = GrpNum
		ATGrpDict[ATGene] = ATNum
		GrpNum += 1
	try: #Look to see if there are any GM orthologs for that gene.
		GMGenesTemp = ATGM[ATGene].split(",")
		for GMGene in GMGenesTemp:
			try: #If there are, see if they already have a group number.
				GMNum = GMGrpDict[GMGene]
				if GMNum != ATNum: #Print an error message if they are in two groups.
					print ("Gene %s is in groups %d and %d.\n" % (GMGene, ATNum, GMNum))
			except KeyError: #If they don't have a group number, give them this number.
				GMGrpDict[GMGene] = ATNum
				#Then look at the AT genes to which the GM gene is orthologous.
				ATListTemp = list(set(GMAT[GMGene].split(",")))
				for ATtempGene in ATListTemp:
					try: #See if those genes have a group number.
						ATGMNum = ATGrpDict[ATtempGene]
						#And print an error message if they are in two groups.
						if ATGMNum != ATNum:
							print ("Gene %s is in groups %d and %d.\n" % (ATtempGene, ATNum, ATGMNum))
					except KeyError: #If they don't, give them this group number
						ATGrpDict[ATtempGene] = ATNum
	except KeyError: #If no GM genes are found for that AT gene, print an error message.
		#This should perhaps be disabled (or sent to a file).
		print ("No GM Genes for %s. \n" % (ATGene))
	#Now do the same for the other four species.
	try:
		OSGenesTemp = ATOS[ATGene].split(",")
		for OSGene in OSGenesTemp:
			try:
				OSNum = OSGrpDict[OSGene]
				if OSNum != ATNum:
					print ("Gene %s is in groups %d and %d.\n" % (OSGene, ATNum, OSNum))
			except KeyError:
				OSGrpDict[OSGene] = ATNum
				ATListTemp = list(set(OSAT[OSGene].split(",")))
				for ATtempGene in ATListTemp:
					try:
						ATOSNum = ATGrpDict[ATtempGene]
						if ATOSNum != ATNum:
							print ("Gene %s is in groups %d and %d.\n" % (ATtempGene, ATNum, ATOSNum))
					except KeyError:
						ATGrpDict[ATtempGene] = ATNum
	except KeyError:
		print ("No OS Genes for %s. \n" % (ATGene))
	try:
		PTGenesTemp = ATPT[ATGene].split(",")
		for PTGene in PTGenesTemp:
			try:
				PTNum = PTGrpDict[PTGene]
				if PTNum != ATNum:
					print ("Gene %s is in groups %d and %d.\n" % (PTGene, ATNum, PTNum))
			except KeyError:
				PTGrpDict[PTGene] = ATNum
				ATListTemp = list(set(PTAT[PTGene].split(",")))
				for ATtempGene in ATListTemp:
					try:
						ATPTNum = ATGrpDict[ATtempGene]
						if ATPTNum != ATNum:
							print ("Gene %s is in groups %d and %d.\n" % (ATtempGene, ATNum, ATPTNum))
					except KeyError:
						ATGrpDict[ATtempGene] = ATNum
	except KeyError:
		print ("No PT Genes for %s. \n" % (ATGene))
	try:
		STGenesTemp = ATST[ATGene].split(",")
		for STGene in STGenesTemp:
			try:
				STNum = STGrpDict[STGene]
				if STNum != ATNum:
					print ("Gene %s is in groups %d and %d.\n" % (STGene, ATNum, STNum))
			except KeyError:
				STGrpDict[STGene] = ATNum
				ATListTemp = list(set(STAT[STGene].split(",")))
				for ATtempGene in ATListTemp:
					try:
						ATSTNum = ATGrpDict[ATtempGene]
						if ATSTNum != ATNum:
							print ("Gene %s is in groups %d and %d.\n" % (ATtempGene, ATNum, ATSTNum))
					except KeyError:
						ATGrpDict[ATtempGene] = ATNum
	except KeyError:
		print ("No ST Genes for %s. \n" % (ATGene))
	try:
		VVGenesTemp = ATVV[ATGene].split(",")
		for VVGene in VVGenesTemp:
			try:
				VVNum = VVGrpDict[VVGene]
				if VVNum != ATNum:
					print ("Gene %s is in groups %d and %d.\n" % (VVGene, ATNum, VVNum))
			except KeyError:
				VVGrpDict[VVGene] = ATNum
				ATListTemp = list(set(VVAT[VVGene].split(",")))
				for ATtempGene in ATListTemp:
					try:
						ATVVNum = ATGrpDict[ATtempGene]
						if ATVVNum != ATNum:
							print ("Gene %s is in groups %d and %d.\n" % (ATtempGene, ATNum, ATVVNum))
					except KeyError:
						ATGrpDict[ATtempGene] = ATNum
	except KeyError:
		print ("No VV Genes for %s. \n" % (ATGene))

for ATGene in ATGrpDict:
	GrpNum = ATGrpDict[ATGene]
	try:
		GrpDict[GrpNum]['AT'] += ATGene+","
	except KeyError:
		GrpDict[GrpNum] = defaultdict(str)
		GrpDict[GrpNum]['AT'] = ATGene+","
for GMGene in GMGrpDict:
	GrpNum = GMGrpDict[GMGene]
	GrpDict[GrpNum]['GM'] += GMGene+","
for OSGene in OSGrpDict:
	GrpNum = OSGrpDict[OSGene]
	GrpDict[GrpNum]['OS'] += OSGene+","
for PTGene in PTGrpDict:
	GrpNum = PTGrpDict[PTGene]
	GrpDict[GrpNum]['PT'] += PTGene+","
for STGene in STGrpDict:
	GrpNum = STGrpDict[STGene]
	GrpDict[GrpNum]['ST'] += STGene+","
for VVGene in VVGrpDict:
	GrpNum = VVGrpDict[VVGene]
	GrpDict[GrpNum]['VV'] += VVGene+","

#Make a list of the sequences in each group, without the final comma in each set of sequences.

GrpDictList = [ ]

for GrpNum in GrpDict.keys():
	OutLine = str(GrpNum)+"\t"+GrpDict[GrpNum]['AT'][:-1]+"\t"+GrpDict[GrpNum]['GM'][:-1]+"\t"+GrpDict[GrpNum]['OS'][:-1]\
	+"\t"+GrpDict[GrpNum]['PT'][:-1]+"\t"+GrpDict[GrpNum]['ST'][:-1]+"\t"+GrpDict[GrpNum]['VV'][:-1]+"\n"
	GrpDictList.append(OutLine)
	
#Write the list to a file.

OutFile = open(OutFileGeneOrth, 'w')
for OutLine in GrpDictList:
	OutFile.write(OutLine)
OutFile.close()
