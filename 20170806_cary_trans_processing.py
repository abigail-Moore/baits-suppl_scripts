#! /usr/bin/env python

#This is a modified version of tbaits_blastn_parse.py and tblast_to_fasta.py that just
#parses the blast output for a single transcriptome and makes files of the full
#sequences and files where each of the sequences is split into various parts.
#This is a modified version of the 20170710 script so that the blast parsing is more stringent.

import sys #to talk to the command line
from collections import defaultdict #to make dictionaries with multiple levels
import numpy #to randomly divide up the sequences

InFileBl = "/home/abby/transcriptomes/combined2_paper1_tree/cary_testing/combined_fasta_out.out"
InFileSeqs = "/home/abby/transcriptomes/combined2_paper1_tree/cary_testing/combined_fasta_in.fa"
LocusFileName = "/home/abby/transcriptomes/general/loci_shortened2.txt"
LocusInfoFileOut = "/home/abby/transcriptomes/combined2_paper1_tree/cary_testing/combined_fasta.info"
SeqFolderOut = "/home/abby/transcriptomes/combined2_paper1_tree/cary_testing/seqs_out/"
SppName = "Cary_Dryad"
NumReps = 10
AlFolder = "/home/abby/transcriptomes/general/c2p1_bbtrees_baitsonly_nn/"
AlFilePre = "c2p1bo_"

#################################################################################

#DictFromFile makes a dictionary from a tab-delimited file, where the keys are in columns KC, and
#the values are in column VC
#from tcontig_selection.py, which was a modified version of the function in from tbaits_intron_removal.py
def DictFromFile(FileName, KC, VC):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[KC]] = Line[VC]
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	#sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict


#OutFileWriting writes an output file from a list of lines to write.
#The lines must already have "\n" at the end.
#from tbaits_intron_removal.py
def OutFileWriting(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line)
	OutFile.close()
	#print("Output file %s written.\n" % (FileName))
	#sys.stderr.write("Output file %s written.\n" % (FileName))

###################################################################################
LocusDict = DictFromFile(LocusFileName, 0, 1) #The dictionary of the locus names from the blast files and the final output files for those names


#looking at the blast output
SeqBlastDict = defaultdict(dict)#SeqBlastDict[SeqName][Locus] = [BitScorelist]
RCDict = {}#RCDict[SeqName] = 'RC'/'no RC'
InFile = open(InFileBl, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	SeqName = Line[0]
	BlHit = Line[1].split("+")[0]
	Locus = LocusDict[BlHit]
	BitScore = float(Line[11])
	if BitScore >= 500:
		try:
			SeqBlastDict[SeqName][Locus].append(BitScore)
		except KeyError:
			SeqBlastDict[SeqName][Locus]=[BitScore]
		#adding the sequence to the list of ones that needs to be reverse-complemented, if it is reversed with respect to the query sequence.
		if int(Line[6]) > int(Line[7]):
			RCDict[SeqName] = "RC"
		else:
			RCDict[SeqName] = "no RC"
InFile.close()

#condensing the lists
SeqInfoDict = defaultdict(dict)
for SeqName in SeqBlastDict.keys():
	for Locus in SeqBlastDict[SeqName].keys():
		SeqInfoDict[SeqName][Locus] = numpy.mean(SeqBlastDict[SeqName][Locus])

#printing the list
SeqNameList = sorted(SeqInfoDict.keys())
OutList = ['Locus\tBlastHits:MeanBitScores\n']
for SeqName in SeqNameList:
	Line = SeqName+"\t"+"\t".join(("%s:%.2f" % (Locus, SeqInfoDict[SeqName][Locus])) for Locus in SeqInfoDict[SeqName].keys())+"\n"
	OutList.append(Line)
OutFileWriting(LocusInfoFileOut,OutList)

#reading the fasta file and making a dictionary
SeqDict = {}
InFile = open(InFileSeqs, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r')
	if Line[0] == ">":
		SeqName = Line[1:]
		SeqDict[SeqName] = ""
	else:
		SeqDict[SeqName] += Line
InFile.close()

#figuring out which sequences we want
SeqsWanted = defaultdict(dict)
OneHit = 0
MultipleHits = 0
for SeqName in SeqInfoDict.keys():
	#If there was only one blast hit for that sequence, go with that locus
	if len(SeqInfoDict[SeqName].keys()) == 1:
		LocusName = SeqInfoDict[SeqName].keys()[0]
		if RCDict[SeqName] == "RC":
			NewSeq2 = Seq(SeqDict[SeqName], IUPAC.unambiguous_dna)
			NewSeq = str(NewSeq2.reverse_complement())
			SeqsWanted[LocusName][SeqName] = NewSeq
		else:
			SeqsWanted[LocusName][SeqName] = SeqDict[SeqName]
		OneHit += 1
	#if there were multiple blast hits, then rank them in order of mean bitscore
	else:
		BestScore = 0.0
		for LocusName in SeqInfoDict[SeqName]:
			if SeqInfoDict[SeqName][LocusName] > BestScore:
				BestLocus = LocusName
				BestScore = SeqInfoDict[SeqName][LocusName]
		if RCDict[SeqName] == "RC":
			NewSeq2 = Seq(SeqDict[SeqName], IUPAC.unambiguous_dna)
			NewSeq = str(NewSeq2.reverse_complement())
			SeqsWanted[LocusName][SeqName] = NewSeq
		else:
			SeqsWanted[LocusName][SeqName] = SeqDict[SeqName]
		MultipleHits += 1
print("Of the %d sequences, %d had one best blast hit and %d had multiple best blast hits.  The best hit according to mean bitscore was chosen for the latter sequences.\n" % (len(SeqInfoDict), OneHit, MultipleHits))

#for LocusName in sorted(SeqsWanted.keys()):
#	print("%s: %s\n" % (LocusName, ",".join(SeqsWanted[LocusName].keys())))

#printing the full sequences
for Locus in SeqsWanted.keys():
	OutFileName = SeqFolderOut+"full_"+SppName+"_"+Locus+".fa"
	OutList = [ ]
	for SeqName in SeqsWanted[Locus]:
		Line = '>'+SeqName+'\n'+SeqsWanted[Locus][SeqName]+'\n'
		OutList.append(Line)
		#print("%s: %s: %d\n" % (Locus, SeqName, len(SeqsWanted[Locus][SeqName])))
	OutFileWriting(OutFileName, OutList)
print("%d full sequence files were written, with names such as %s.\n" % (len(SeqsWanted.keys()), OutFileName))

#making a bash script to align them
OutFileName = SeqFolderOut+"alignment_script.sh"
OutList = ["#! /bin/bash\n"]
for Locus in SeqsWanted.keys():
	Line = "mafft --localpair --add "+SeqFolderOut+"full_"+SppName+"_"+Locus+".fa --thread -1 "+AlFolder+AlFilePre+Locus+".fa > "+SeqFolderOut+"full_"+SppName+"_"+Locus+"_al.fa\n"
	OutList.append(Line)
OutFileWriting(OutFileName, OutList)
print("The script %s must now be executed to align the sequences to the backbone alignments.\n" % (OutFileName))