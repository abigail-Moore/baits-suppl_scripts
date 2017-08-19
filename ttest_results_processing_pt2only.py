#! /usr/bin/env python

#three files to parse:
#SeqInfo_Group1.txt with FileName [tab] ContigName [tab] ContigLength--produced by tseq_lengths.py
#pt1cb_Contig_Fates.txt with Locus [tab] Group [tab] Contig [tab] Contig_Fate (and a header row) in the pt1_combined folder.
#where Contig_Fate is the combined sequence name the contig was put into (prefaced with good_) or, if the contig could not be classified,
#either undiv_, meaning it was used, or not_using_, meaning it was not.
#The Group in this case is all the same (all Lophophora williamsii), so it can be ignored.
#pt1gt1_combined_seqs.txt with Locus [tab] IndName [tab] New_ParalogName [tab] Paralogs_Combined (comma separated) (and  header row) in the pt1_genefams1 folder.

import sys #to talk to the command line
from collections import defaultdict #dictionaries with multiple levels

Usage = '''
ttest_results_processing.py version 1.0 16 July 2017
This script assesses the accuracy of contig classification by looking at the 
Contig_Fates file and also the lengths of the various contigs.
ttest_results_processing.py
[name of the file with sequence name [tab] sequence length]
[name of the Contig_Fates.txt file]
[cutoff value for a sequence to be considered reliably classifiable]
[name of the combined_seqs.txt file]
[file name for the list of individuals--expecting the individual name in the second column, as if it is a list of individuals in groups]
[prefix for the output files, including the folder]
[the file name for the cummulative stats file]
'''

print("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 8:
	sys.exit("ERROR!  This script requires 7 additional arguments and you supplied %d.\n %s" % (len(sys.argv)-1, Usage))
ContigLenFN = sys.argv[1]
ContigFateFN = sys.argv[2]
CutoffValue = int(sys.argv[3])
CombSeqFN = sys.argv[4]
IndListFN = sys.argv[5]
OutFilePre = sys.argv[6]
StatsFN = sys.argv[7]

#################################################################

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

#ListDictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value.  Each key has multiple values, and they
#will be made into a list
#modified from tparalaog_combiner.py
def ListDictFromFile(FileName, KeyCol, ListCol):
	TempDict = defaultdict(list)
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[KeyCol]].append(Line[ListCol])
	InFile.close()
	for Key in TempDict:
		ListTemp = TempDict[Key]
		ListTemp = list(set(ListTemp))
		TempDict[Key] = ListTemp
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
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
	print("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	#sys.stderr.write("Column number %d was read from the file %s.\n" % (ColNum+1, FileName))
	return TempList
	#This is IndList

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

######################################################################

#reading the input files
ContigLenDict = DictFromFile(ContigLenFN, 0, 1) #ContigLenDict[ContigName] = ContigLen
IndList = CaptureColumn(IndListFN, 1)
ContigFateDict = defaultdict(dict)#ContigFateDict[Locus][ContigName] = ContigFate
CombSeqDict = defaultdict(dict) #CombSeqDict[Locus][NewPN--New Paralog Name] = [CombPars]--list of paralogs that were combined
SeqCombiningDict = defaultdict(dict) #SeqCombiningDict[Locus][PN] = NewPN
#SeqCombiningDictFinal = defaultdict(dict) #SeqCombiningDictFinal[Locus][PN] = NewPN--this one has the original paralog name with the final new paralog name
#S3ContigFateDict = defaultdict(dict)#S3ContigFateDict[Locus][ContigName] = ContigFate (ContigFate after being combined in step3)

#Filling out ContigFateDict from ContigFateFN
InFile = open(ContigFateFN, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	if Line[0] != "Locus":
		Locus = Line[0]
		ContigName = Line[2]
		ContigFate = Line[3]
		if "_to_" not in ContigName: #We can't deal with split fragments, because they won't be in the list of contig lengths.
			ContigFateDict[Locus][ContigName] = ContigFate
InFile.close()

'''
#filling out CombSeqDict from CombSeqFN ([pre]gt1_combined_seqs.txt)
#pt1gt1_combined_seqs.txt with Locus [tab] IndName [tab] New_ParalogName [tab] Paralogs_Combined (comma separated) (and header row).
InFile = open(CombSeqFN, 'rU')
for Line in InFile:
	Line = Line.strip('\n').strip('\r').split('\t')
	if Line[0] != "Locus":
		Locus = Line[0]
		IndName = Line[1]
		NewPN = Line[2].split('.')[1]
		CombPars = [Item.split('.')[1] for Item in Line[3].split(',')]
		if IndName in IndList:
			CombSeqDict[Locus][NewPN] = CombPars
			for PN in CombPars:
				SeqCombiningDict[Locus][PN] = NewPN
InFile.close()

for Locus in ContigFateDict:
	for ContigName in ContigFateDict[Locus]:
		PN = ContigFateDict[Locus][ContigName]
		#Reformatting the paralog name to follow the SeqCombiningDict
		PNTemp = ("_".join(PN.split("_")[1:]))
		try:
			S3ContigFateDict[Locus][ContigName] = SeqCombiningDict[Locus][PNTemp]
			#print("New paralog %s for contig %s, instead of %s.\n" % (SeqCombiningDict[Locus][PNTemp], ContigName, PNTemp))
		except KeyError:
			S3ContigFateDict[Locus][ContigName] = PN
'''

#figuring out how to parse the files
TotalContigs = 0
SeqFateDict = defaultdict(dict) #SeqFateDict[Locus][SeqName][['ClassifList'] = [list of ASName]['NumSeqs']
AssembledSeqDict = defaultdict(dict) #AssembeledSeqDict[Locus][ASName][['ContigList']['NumContigs']]
SeqList = [ ]
for Locus in ContigFateDict:
	for ContigName in ContigFateDict[Locus]:
		if int(ContigLenDict[ContigName]) > CutoffValue:
			SeqName = (ContigName.split("-")[1]).split(".")[0]
			SeqList.append(SeqName)
			ASName = ContigFateDict[Locus][ContigName]
			TotalContigs += 1
			try:
				SeqFateDict[Locus][SeqName]['ClassifListAll'].append(ASName)
				SeqFateDict[Locus][SeqName]['NumSeqsAll'] += 1
			except KeyError:
				SeqFateDict[Locus][SeqName] = {}
				SeqFateDict[Locus][SeqName]['ClassifListAll'] = [ASName]
				SeqFateDict[Locus][SeqName]['ClassifListPruned'] = [ ]
				SeqFateDict[Locus][SeqName]['NumSeqsAll'] = 1
				SeqFateDict[Locus][SeqName]['NumSeqsPruned'] = 0
			#print("%s: %s" % (ASName, ASName.split("_")[0]))
			if ASName.split("_")[0] not in ["not", "undiv"]:
				SeqFateDict[Locus][SeqName]['ClassifListPruned'].append(ASName)
				SeqFateDict[Locus][SeqName]['NumSeqsPruned'] += 1
				try:
					AssembledSeqDict[Locus][ASName]['ContigList'].append(SeqName)
					AssembledSeqDict[Locus][ASName]['NumContigs'] += 1
				except KeyError:
					AssembledSeqDict[Locus][ASName] = { }
					AssembledSeqDict[Locus][ASName]['ContigList'] = [SeqName]
					AssembledSeqDict[Locus][ASName]['NumContigs'] = 1
SeqList = list(set(SeqList))


#condensing the lists for SeqFateDict
GoodClassifPruned = 0
GoodClassifContigsPruned = 0
BadClassifPruned = 0
BadClassifContigsPruned = 0
GoodClassifAll = 0
GoodClassifContigsAll = 0
BadClassifAll = 0
BadClassifContigsAll = 0
for Locus in SeqFateDict:
	for SeqName in SeqFateDict[Locus]:
		SeqFateDict[Locus][SeqName]['ClassifListPruned'] = list(set(SeqFateDict[Locus][SeqName]['ClassifListPruned']))
		SeqFateDict[Locus][SeqName]['ClassifListAll'] = list(set(SeqFateDict[Locus][SeqName]['ClassifListAll']))
		if len(SeqFateDict[Locus][SeqName]['ClassifListPruned']) == 1:
			GoodClassifPruned += 1
			GoodClassifContigsPruned += SeqFateDict[Locus][SeqName]['NumSeqsPruned']
		else:
			BadClassifPruned += 1
			BadClassifContigsPruned += SeqFateDict[Locus][SeqName]['NumSeqsPruned']
		if len(SeqFateDict[Locus][SeqName]['ClassifListAll']) == 1:
			GoodClassifAll += 1
			GoodClassifContigsAll += SeqFateDict[Locus][SeqName]['NumSeqsAll']
		else:
			BadClassifAll += 1
			BadClassifContigsAll += SeqFateDict[Locus][SeqName]['NumSeqsAll']
print("For all %d fragments that are longer than %d bases, %d were assembled correctly (from %d original sequences), while %d were assembled incorrectly (from %d original sequences), and %d were not assembled.\n" % (TotalContigs, CutoffValue, GoodClassifContigsPruned, GoodClassifPruned, BadClassifContigsPruned, BadClassifPruned, TotalContigs-GoodClassifContigsPruned-BadClassifContigsPruned))
print("For all %d fragments that are longer than %d bases, %d were assembled correctly (from %d original sequences), while %d were assembled incorrectly (from %d original sequences), when considering all fragments.\n" % (TotalContigs, CutoffValue, GoodClassifContigsAll, GoodClassifAll, BadClassifContigsAll, BadClassifAll))
StatsFile = open(StatsFN, 'a')
StatsFile.write(("\t").join([ContigFateFN, str(CutoffValue), str(TotalContigs), str(GoodClassifContigsPruned), str(BadClassifContigsPruned), str(TotalContigs-GoodClassifContigsPruned-BadClassifContigsPruned), str(len(SeqList)), str(GoodClassifPruned), str(BadClassifPruned), str(len(SeqList)-GoodClassifPruned-BadClassifPruned)])+"\t")


#condensing the lists for AssembledSeqDict
GoodClassif = 0
GoodClassifSeqs = 0
BadClassif = 0
BadClassifSeqs = 0
for Locus in AssembledSeqDict:
	for ASName in AssembledSeqDict[Locus]:
		AssembledSeqDict[Locus][ASName]['ContigList'] = list(set(AssembledSeqDict[Locus][ASName]['ContigList']))
		if len(AssembledSeqDict[Locus][ASName]['ContigList']) == 1:
			GoodClassif += 1
			GoodClassifSeqs += AssembledSeqDict[Locus][ASName]['NumContigs']
		else:
			BadClassif += 1
			BadClassifSeqs += AssembledSeqDict[Locus][ASName]['NumContigs']
print("Of the %d assembled sequences, %d were composed of fragments from the same original sequence (from %d of the original sequences), while %d were composed of fragments from different original sequences (from %d original sequences).\n" % (GoodClassifSeqs+BadClassifSeqs, GoodClassifSeqs, GoodClassif, BadClassifSeqs, BadClassif))
StatsFile.write(("\t").join([str(GoodClassifSeqs+BadClassifSeqs), str(GoodClassifSeqs), str(BadClassifSeqs)])+"\n")
StatsFile.close()


#writing the three output files
OutFileName = OutFilePre+"Seq_Fate_All_Contigs_"+str(CutoffValue)+".txt"
OutList = ['Locus\tName_of_Original_Sequence\tNumber_of_Assembled_Sequences\tList_of_Assembled_Sequences\n']
for Locus in sorted(SeqFateDict.keys()):
	for SeqName in sorted(SeqFateDict[Locus].keys()):
		Line = ("\t").join([Locus, SeqName, str(SeqFateDict[Locus][SeqName]['NumSeqsAll']), (",").join(SeqFateDict[Locus][SeqName]['ClassifListAll'])])+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)

OutFileName = OutFilePre+"Seq_Fate_Pruned_Contigs_"+str(CutoffValue)+".txt"
OutList = ['Locus\tName_of_Original_Sequence\tNumber_of_Assembled_Sequences\tList_of_Assembled_Sequences\n']
for Locus in sorted(SeqFateDict.keys()):
	for SeqName in sorted(SeqFateDict[Locus].keys()):
		Line = ("\t").join([Locus, SeqName, str(SeqFateDict[Locus][SeqName]['NumSeqsPruned']), (",").join(SeqFateDict[Locus][SeqName]['ClassifListPruned'])])+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)

OutFileName = OutFilePre+"Assembled_Sequence_Information_"+str(CutoffValue)+".txt"
OutList = ['Locus\tAssembled_Sequence_Name\tNumber_of_Contigs_in_Sequence\tList_of_Template_Sequences\n']
for Locus in sorted(AssembledSeqDict.keys()):
	for ASName in sorted(AssembledSeqDict[Locus].keys()):
		Line = ("\t").join([Locus, ASName, str(AssembledSeqDict[Locus][ASName]['NumContigs']), (",").join(AssembledSeqDict[Locus][ASName]['ContigList'])])+"\n"
		OutList.append(Line)
OutFileWriting(OutFileName, OutList)