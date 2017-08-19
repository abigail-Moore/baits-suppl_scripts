
##########################################
Methods summary for designing probes from transcriptomes and scripts used
4 June 2014

##Procedure for loci from the transcriptome:

1.  Build clusters of orthologous genes from six model plant species (Arabidopsis thaliana, Glycine max, Oryza sativa, Populus trichocarpa, Solanum tuberosum, Vitis vinifera) using the list of groups of orthologous genes from Ensemble.  These clusters are non-overlapping (i.e., no gene belongs to multiple clusters).
(script: gene_ortho_2.py to make the lists of genes and ref_groups_2.py to make the fasta files for each cluster)

2.  Determine which GO terms apply to each cluster using the GO terms given to the Arabidopsis genes in the cluster and classify each GO term as belonging to all or only a subset of Arabidopsis genes that make up the cluster.
(script: ortho_GO_1.py)

3.  Use BLAST to align our transcriptome sequences, the genes from the Beta vulgaris genome, and the 1KP Mollugo sequences to the clusters.

4.  Go through the BLAST results, determine which sequences belong to each cluster, and add the sequences to the end of the fasta file for that cluster.
(script: blast_dict_2.py)

5.  Determine whether each cluster has at least one Mollugo sequence, at least one cactus sequence, and at least one sequence from another member of the Portulacineae and throw out the remaining clusters.
(script: tgroup_filter_2osc.py)

6.  Use Muscle (or in rare cases when Muscle was taking a prohibitively long time Clustal) to align the genes in each cluster (on Oscar).

7.  Look through each alignment.  If the average sequence divergence is less than 40%, the cluster is ready for final analysis (step 8).  If the average sequence divergence is more then 40%, make a tree, root the tree at the midpoint, and examine each subtree.  If the subtree has at least one Mollugo sequence, at least one cactus sequence, and at least one sequence from another member of the Portulacineae, make a new sequence file and go back to step 6 to align it again.  If not, throw it out.
(scripts: tree_parse.r, tree_parse_foros.py)

8.  Look through the finished alignments that have fulfilled all of the above criteria and assign GO terms to the alignments that have Arabidopsis loci.  (We actually did not end up using any of the alignments with no Arabidopsis loci.)
(script: talign_groups.py)

9.  Build trees for each alignment with RAxML so we have a better idea how the sequences are related to each other.

9.  Based on the GO terms of each alignment, select alignments of potential interest for designing baits.
(potential scripts: tGOterm_list_parse_2.py—makes a list of the alignments with certain GO terms, tseq_length.py—makes a list of all of the sequences in a set of alignments with their lengths, the alignments of interest can be marked on this list)

10.  Choose the specific sequences in the alignments from which baits should be designed.  Alignments were chosen according to coverage and phylogenetic position (trying to span the entire clade).

11.  Make new a new fasta file for each locus that just includes the sequences from which we want to design baits.  Align these sequences using muscle.
(script: talbaits_2.py)

##Procedure for the loci for which we already had alignments (from P.-A. Christin, M. Arakaki, C. P. Osborne, and E. J. Edwards.  2015.  Genetic enablers underlying the clustered evolutionary origins of C4 photosynthesis in angiosperms.  Molecular Biology and Evolution 32:846-858.):

1.  Combine all alignments into a gigantic file for BLAST searching.
(script: talignment_combine.py)

2.  BLAST the 1KP Mollugo transcriptomes against the alignments.

3.  Add the Mollugo sequences to the ends of the alignments to which they BLAST.
(script: blast_dict_2.py)

4.  Align the sequences using Muscle and build trees using RAxML.

5.  Look through the trees to separate the paralogues and remove very divergent Mollugo sequences.
(scripts: tseq_length.py—makes the file to be annotated with the fates of the sequences, talsplit.py—splits the alignments according to the annotations in the file)

6.  Align the new sets of sequences using Muscle and build trees using RAxML.

7.  Look through the smaller alignments and choose sequences for bait design.

8.  Make new a new fasta file for each locus that just includes the sequences from which we want to design baits.  Align these sequences using muscle.
(script: talbaits_2.py)