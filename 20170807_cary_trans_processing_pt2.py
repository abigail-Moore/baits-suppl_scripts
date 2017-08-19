#! /usr/bin/env python

#This is a modified version of tbaits_blastn_parse.py and tblast_to_fasta.py that just
#parses the blast output for a single transcriptome and makes files of the full
#sequences and files where each of the sequences is split into various parts.
#This is a modified version of the 20170710 script so that the blast parsing is more stringent.

import sys #to talk to the command line
from collections import defaultdict #to make dictionaries with multiple levels
import numpy #to randomly divide up the sequences
import subprocess #We want to be able to talk to the command line.

LocusListFN = "/home/abby/transcriptomes/general/Locus_List_all2.txt"
SeqFolderIn = "/home/abby/transcriptomes/combined2_paper1_tree/cary_testing/seqs_out_nn/"
SeqFilePre = "full_Cary_Dryad_"
SeqFolderOut = "/home/abby/transcriptomes/combined2_paper1_tree/cary_testing/seqs_out_parts/"
OutFilePre = "Cary"
GIListFN = "/home/abby/transcriptomes/combined2_paper1_tree/cary_testing/cary_Ind_List_Groups.txt"
GListFN = "/home/abby/transcriptomes/combined2_paper1_tree/cary_testing/cary_Group_List.txt"
NumReps = 10

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

#CaptureColumn makes a list from a specified column in a file.  This is useful
#for reusing various text files that previous programs needed.
#The column number has to follow python numbering, so 0 for the first, 1 for the
#second, etc.
#from tseq_placer_dup.py
def CaptureColumn(FileName, ColNum):
	TempList = [ ]
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\n').strip('\r').split('\t')
		TempList.append(Line[ColNum])
	InFile.close()
	#print("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	#sys.stderr.write("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	return TempList
	#This is LocusList


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

GroupDict = DictFromFile(GListFN, 0, 1) #GroupDict[Group] = Abbr
GroupIndDict = DictFromFile(GIListFN, 1, 0) #GroupIndDict[IndName] = Group
LocusList = CaptureColumn(LocusListFN, 0)

#making the new folders
OutLine = "mkdir "+SeqFolderOut
subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
for Group in GroupDict:
	OutLine = "mkdir "+SeqFolderOut+Group
	subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
	OutLine = "mkdir "+SeqFolderOut+Group+"/s2_final_"+GroupDict[Group]
	subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]
	for SetNum in range(1,NumReps+1):
		OutLine = "mkdir "+SeqFolderOut+Group+"/s2_final_"+GroupDict[Group]+"/cary_grp"+str(SetNum)+"_exons"
		subprocess.Popen(OutLine, shell=True, stdout=subprocess.PIPE).communicate()[0]

SeqDict = defaultdict(dict) #SeqDict[Locus][SeqName] = Seq
LocusLenDict = { }
for Locus in LocusList:
	InFileName = SeqFolderIn+SeqFilePre+Locus+"_al.fa"
	try: 
		InFile = open(InFileName, 'rU')
		for Line in InFile:
			Line = Line.strip('\n').strip('\r')
			if Line[0] == ">":
				SeqName = Line[1:]
				SeqDict[Locus][SeqName] = ""
			else:
				SeqDict[Locus][SeqName] += Line
				LocusLenDict[Locus] = len(SeqDict[Locus][SeqName])
		InFile.close()
	except IOError:
		print("Locus %s does not have a sequence file.\n" % (Locus))

GroupInfoDict = defaultdict(dict) #GroupInfoDict[SetNum][FragName] = SeqLen
for Locus in LocusLenDict:
	SeqLen = LocusLenDict[Locus]
	NumBreaks = SeqLen/300
	#print("\n\n%s: %d: %d" % (Locus, SeqLen, NumBreaks))
	if NumBreaks > 0:
		for SetNum in range(0,NumReps):
			OutSeqDict = defaultdict(dict) #OutSeqDict[Group][FragName] = SeqLen
			BreakPoints = numpy.random.choice(SeqLen, NumBreaks, replace=False)
			BreakPoints = sorted(BreakPoints)
			NumSeqs = 0
			for SeqName in SeqDict[Locus]:
				IndName = SeqName.split("-")[0]
				Group = GroupIndDict[IndName]
				NumSeqs += 1
				#print("%s: %s: %s" % (SeqName, IndName, Group))
				Seq = SeqDict[Locus][SeqName]
				for Break in range(0,NumBreaks):
					FragName = SeqName+"."+str(Break+1)
					if Break == 0:
						NewFrag = Seq[0:BreakPoints[Break]]
					else:
						NewFrag = Seq[BreakPoints[Break-1]:BreakPoints[Break]]
					NewFrag = NewFrag.replace("-","")
					FragLen = len(NewFrag)
					if FragLen > 1:
						OutSeqDict[Group][FragName] = NewFrag
						GroupInfoDict[SetNum][FragName] = FragLen
				#print("%s: %s: %s" % (SeqName, IndName, Group))
			#now to print the sequences
			for Group in OutSeqDict:
				#print("%s: %d" % (Group, SetNum+1))
				OutFileName = SeqFolderOut+Group+"/s2_final_"+GroupDict[Group]+"/cary_grp"+str(SetNum+1)+"_exons"+"/grp"+str(SetNum+1)+"e1_"+Locus+".fa"
				OutList = [ ]
				for FragName in OutSeqDict[Group]:
					Line = ">"+FragName+"\n"+OutSeqDict[Group][FragName]+"\n"
					OutList.append(Line)
				OutFileWriting(OutFileName, OutList)
print("%d sets of files were made, for %d loci.\n" % (NumReps, len(LocusLenDict)))


#print the GroupInfoDict for each group
for SetNum in range(0,NumReps):
	OutFileName =  SeqFolderOut+"SeqInfo_Group"+str(SetNum+1)+".txt"
	OutList = [ ]
	for SeqName in sorted(GroupInfoDict[SetNum].keys()):
		Line = SeqName+"\t"+str(GroupInfoDict[SetNum][SeqName])+"\n"
		OutList.append(Line)
	OutFileWriting(OutFileName, OutList)
print("%d files of information about the fragmented sequences, with names such as %s, were written.\n" % (NumReps, OutFileName))