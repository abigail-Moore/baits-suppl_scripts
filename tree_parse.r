#! /usr/bin/Rscript

#This runs rather quickly, so I just run it and the subsequent python script on my
#laptop, and only run the alignments (which take a long time) on Oscar.

rm(list=ls())

#We will need ape to do phylogenetically-related calculations.
library(ape)
library(phangorn)

#tree_parse.r version 1.0 4 April 2014 Abby Moore
#This script looks through a set of alignments and calculates the distance matrix
#between the sequences.  If the mean distance is below a given threshold then
#it considers the tree to be ready for analysis.  If the mean distance is above
#that threshold, then it splits the tree at the root and looks at the two subtrees.
#The subtrees that have at least one Mollugo, cactus, and other Portulacineae are
#put in a list for further analysis and the list of sequences is saved to a file.
#This script is called as follows:
#tree_parse.r [Working Directory] [InFile] [OutFile]

args <- commandArgs(TRUE)
WorkingDirectory <- args[1]
InFile <- args[2]
OutFile <- args[3]

#setwd("/home/abby/transcriptomes/newgroups")
setwd(WorkingDirectory)

#Reading the names of the alignment files from the input file.
#We do not want the file names to be read as factors (thus as.is = TRUE).
#SeqList <- read.table("nGroup_RList.txt", header=F, as.is = TRUE)
SeqList <- read.table(InFile, header=F, as.is = TRUE)

for (i in 1:nrow(SeqList)) {
	FileName <- SeqList[i,1]
	#print(FileName)
	GroupNum <- SeqList[i,2]
	#print(GroupNum)
	#Calculating the distance matrix.
	Align1 <- read.dna(FileName, format = "fasta")
	DMatrix <- dist.dna(Align1, model = "K80", pairwise.deletion = TRUE)
	DMean <- mean(DMatrix)
	DMax <- max(DMatrix)
	#In order to be able to calculate the tree, all of the sequences need to
	#have defined distances between them.  For that reason, I substitute the
	#distances NaN and Inf with the maximum distance between sequences.  I do
	#not know if this is a good idea, but I hope it will at least allow the
	#sequences to fall onto the correct halves of the tree, and thus enable
	#subsequent, better alignments
	if (DMean == "NaN" | DMean == "Inf"){
		DMax2 <- max(DMatrix, na.rm = TRUE)
		if (DMax2 != "Inf") {
			for (j in 1:length(DMatrix)){
				if (DMatrix[j] == "NaN") {
					DMatrix[j] <- DMax2
				}
			}
		}
		else {
			for (j in 1:length(DMatrix)){
				if (DMatrix[j] == "Inf") {
					DMatrix[j] <- NaN
				}
			}
			DMax2 <- max(DMatrix, na.rm = TRUE)
			for (j in 1:length(DMatrix)) {
				if (DMatrix[j] == "NaN") {
					DMatrix[j] <- DMax2
				}
			}
		}
		DMean <- mean(DMatrix)
	}
	#Here I am using the arbitrary mean distance threshold of 0.4 to decide whether
	#to keep a tree as is for analysis or to split it further.	
	if (DMean > 0.40) {
		SeqList[i,3] <- "split"
		#making the full tree (tree1, just saved as nGroup_GroupNum.tre)
		tree1 <- upgma(DMatrix)
		midpoint(tree1)
		tFileName <- paste("nGroup_", GroupNum, ".tre", sep = "", collapse = NULL)
		write.tree(tree1, tFileName)
		NumTaxa = length(tree1$tip.label)
		#making and examining the first subtree (tree2)
		tree2 <- extract.clade(tree1, node = NumTaxa+2)
		tFileName <- paste("nGroup_", GroupNum, "_A.tre", sep = "", collapse = NULL)
		write.tree(tree2, tFileName)
		#Making sure we have at least one cactus, at least one Mollugo, and at least one
		#other Portulacineae in this tree.
		NumCact <- 0
		NumMoll <- 0
		NumPort <- 0
		for (Taxon in tree2$tip.label){
			TipLabel <- strsplit(Taxon, "_")
			if (TipLabel[[1]][1] == "C") NumCact <- NumCact + 1
			if (TipLabel[[1]][1] == "M") NumMoll <- NumMoll + 1
			if (TipLabel[[1]][1] == "P") NumPort <- NumPort + 1
		}
		if ((NumCact != 0) && (NumMoll != 0) && (NumPort != 0)) {
			SeqList[i,4] <- "continue"
			tFileName <- paste("nGroup_", GroupNum, "_Alist.txt", sep = "", collapse = NULL)
			write(tree2$tip.label, tFileName)
		}
		else {
			Message <- paste ("rejected, NumCact: ", NumCact, " NumMoll: ", NumMoll, " NumPort: ", NumPort, sep = "", collapse = NULL)
			SeqList[i,4] <- Message
		}
		#making and examining the second subtree (tree3)
		for (j in 1:dim(tree1$edge)[1]) {
			if (tree1$edge[j,1] == NumTaxa + 1)
			RootTree3 <- tree1$edge[j,2]
		}
		tree3 <- extract.clade(tree1, node = RootTree3)
		#If midpoint rooting resulted in one sequence being sister to all the rest,
		#the two subtrees will be the same.  In this case, we are only keeping one of them.
		if ((length(tree3$tip.label) == length(tree2$tip.label)) && (length(tree3$tip.label) != 0.5*length(tree1$tip.label))) {
			SeqList[i,5] <- "same as tree A"
		}
		#If the two subtrees are different, doing the same thing for tree B as we did for tree A.
		else {
			tFileName <- paste("nGroup_", GroupNum, "_B.tre", sep = "", collapse = NULL)
			write.tree(tree3, tFileName)
			NumCact <- 0
			NumMoll <- 0
			NumPort <- 0
			for (Taxon in tree3$tip.label){
				TipLabel <- strsplit(Taxon, "_")
				if (TipLabel[[1]][1] == "C") NumCact <- NumCact + 1
				if (TipLabel[[1]][1] == "M") NumMoll <- NumMoll + 1
				if (TipLabel[[1]][1] == "P") NumPort <- NumPort + 1
			}
			if ((NumCact != 0) && (NumMoll != 0) && (NumPort != 0)) {
				SeqList[i,5] <- "continue"
				tFileName <- paste("nGroup_", GroupNum, "_Blist.txt", sep = "", collapse = NULL)
				write(tree3$tip.label, tFileName)
			}
			else {
				Message <- paste ("rejected, NumCact: ", NumCact, " NumMoll: ", NumMoll, " NumPort: ", NumPort, sep = "", collapse = NULL)
				SeqList[i,5] <- Message
			}
		}
	}
	#If the mean distance between sequences is less than 0.40, the tree is ready for analysis.
	else {
		SeqList[i,3] <- "ready"
	}
}

#All of this information needs to be printed to a file.
write.table(SeqList, OutFile, quote = FALSE, sep = "\t", row.names = FALSE)
