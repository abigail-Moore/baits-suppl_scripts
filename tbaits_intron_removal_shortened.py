#! /usr/bin/env python

#tbaits_intron_removal.py version 1.1 17 Nov. 2015 Abby Moore
#This is a shortened version of tbaits_intron_removal.py that just makes the output script, given a list of loci,
#for use when introns do not need to be removed.
#version 1.1 has been changed to take LocusList as input and not rely on the BlastFileList
#File formats:
#BlastFileList:
'''BlastFileName [0] (bf2_nhd.out)
SeqFileName [1] (sf2_nhd.fa)'''
#LocusKeyFile:
'''PLocus [0] (apl1)
GLocus [1] (apl)'''


#examples:
'''
tbaits_intron_removal.py /home/abby/transcriptomes/TS31_1/sandbox/gfams/loci_shortened_ppc.txt bf2_ sf2_ ~/transcriptomes/TS31_1/sandbox/gfams/gfcontigs/ same ~/transcriptomes/TS31_1/sandbox/gfams/gfouttemp/ of3t_ ~/transcriptomes/baits_Bv/Bv_groups/outgroup_list.txt ~/transcriptomes/baits_Bv/Bv_groups/ none none none
tbaits_intron_removal.py LocusKeyFile BlastFilePre SeqFilePre SeqFolder BlastFolder OutFolder OutFilePre OGFileName AlFolder AlFilePre AlFilePost FilePath
tbaits_intron_removal_shortened.py LocusListFileName OutFolder OutFilePre OGFileName AlFolder AlFilePre AlFilePost FilePath
'''

from collections import defaultdict
import sys

Usage = '''
tbaits_intron_removal_shortened.py version 1.0
This is a shortened version of tbaits_intron_removal.py that just makes the
final script, without removing introns.
tbaits_intron_removal_shortened.py
[file name for list of loci]
[folder where output files should be written, or "same", if it is the same 
folder as where the sequences are found]
[prefix for output files, or "none", if none]
[for phylogenetic analysis: file with the list of outgroups for each locus]
[folder containing the template alignments and backbone trees]
[alignment file prefix or "none"]
[alignment file ending or "none"]
[path to the fasta_to_phylip.py script or "none" if it is in the default path]
'''

print ("%s\n" % (" ".join(sys.argv)))

if len(sys.argv) != 9:
	sys.exit("ERROR!  This script requires 8 additional arguments, and you supplied %d.  %s\n" % (len(sys.argv)-1, Usage))
LocusListFileName = sys.argv[1]
OutFolder = sys.argv[2]
if OutFolder == "same":
	OutFolder = SeqFolder
OutFilePre = sys.argv[3]
if OutFilePre == "none":
	OutFilePre = ""
OGFileName = sys.argv[4]
AlFolder = sys.argv[5]
AlFilePre = sys.argv[6]
if AlFilePre == "none":
	AlFilePre = ""
AlFilePost = sys.argv[7]
if AlFilePost == "none":
	AlFilePost = ""
FilePath = sys.argv[8]
if FilePath[-1] != "/":
	FilePath += "/"
if FilePath == "none/":
	FilePath = ""


#adding slashes, if necessary:
if OutFolder[-1] != "/":
	OutFolder += "/"

###########functions

#CaptureColumn makes a list from a specified column in a file.  This is useful
#for reusing various text files that previous programs needed.
#The column number has to follow python numbering, so 0 for the first, 1 for the
#second, etc.
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
	#This is LocusList

#DictFromFile makes a dictionary from a tab-delimited file, where the first
#column is the key and the second is the value
def DictFromFile(FileName):
	TempDict = { }
	InFile = open(FileName, 'rU')
	for Line in InFile:
		Line = Line.strip('\r').strip('\n').split('\t')
		TempDict[Line[0]] = Line[1]
	InFile.close()
	print("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	sys.stderr.write("%d lines were read from the file %s and saved to a dictionary.\nExample: %s: %s\n" % (len(TempDict), FileName, Line[0], Line[1]))
	return TempDict
	#This is OutGroupDict

#OutFileWriting writes an output file from a list of lines to write.
#The lines must already have "\n" at the end.
def OutFileWriting(FileName, MyList):
	OutFile = open(FileName, 'w')
	for Line in MyList:
		OutFile.write(Line)
	OutFile.close()
	print("Output file %s written.\n" % (FileName))
	sys.stderr.write("Output file %s written.\n" % (FileName))
		
#MRScriptWriter writes an output script for analysis of the sequences using mafft and raxml
def MRScriptWriter(SeqFileDict, Folder, Prefix, OGDict, AFolder, APre, APost, Path):
	for LocusGroup in SeqFileDict:
		#We need two OutLists because otherwise the files can be removed as the new files are being created.
		OutList = [ ]
		OutList2 = [ ]
		OutLocusList = [ ]
		if LocusGroup != "Ambig":
			for Locus in SeqFileDict[LocusGroup]:
				NamePart = SeqFileDict[LocusGroup][Locus]
				Line = "rm "+Folder+NamePart+"_exons_al.fa && "
				Line += "rm "+Folder+"RAxML*"+NamePart+"\n"
				OutList.append(Line)
				Line = "mafft --addfragments "+Folder+NamePart+".fa --quiet --thread -1 "+AFolder+APre+Locus+APost+".fa > "+Folder+NamePart+"_exons_al.fa && "
				Line += Path+"fasta_to_phylip.py "+Folder+NamePart+"_exons_al.fa && "
				Line += "raxmlHPC -f v -s "+Folder+NamePart+"_exons_al.phy -n "+NamePart+" -t "+AFolder+"RAxML_bipartitions."+APre+Locus+" -m GTRCAT -o "+OGDict[Locus]+" -w "+Folder+"\n"
				OutList2.append(Line)
				OutLocusList.append(Locus+"\n")
			OutFileName = Folder+Prefix+LocusGroup+"Analysis_Script.sh"
			OutList += OutList2
			OutFileWriting(OutFileName, OutList)
			print("The shell script for analyzing this %s group of sequences data was written to %s.\n" % (LocusGroup, OutFileName))
			sys.stderr.write("The shell script for analyzing this %s group of sequences data was written to %s.\n" % (LocusGroup, OutFileName))
			OutFileName = Folder+Prefix+LocusGroup+"Locus_List.txt"
			OutFileWriting(OutFileName,OutLocusList)
			print("The list of loci was written to the file %s.\n" % (OutFileName))
			sys.stderr.write("The list of loci was written to the file %s.\n" % (OutFileName))
		else:
			for Locus in SeqFileDict[LocusGroup]:
				Line = Locus+"\t"+SeqFileDict[LocusGroup][Locus]+"\n"
				OutList.append(Line)
			OutFileName = Folder+Prefix+LocusGroup+"_Ambiguous_Contigs.txt"
			OutFileWriting(OutFileName, OutList)
			print("The list of files of ambiguous sequences was written to %s.\n" % (OutFileName))
			sys.stderr.write("The list of files of ambiguous sequences was written to %s.\n" % (OutFileName))

##########################################################################################################################3


LocusList = CaptureColumn(LocusListFileName, 0)
OutGroupDict = DictFromFile(OGFileName)

FileDict = defaultdict(dict)
for Locus in LocusList:
	FileDict[''][Locus] = OutFilePre+Locus
MRScriptWriter(FileDict, OutFolder, OutFilePre, OutGroupDict, AlFolder, AlFilePre, AlFilePost, FilePath)
