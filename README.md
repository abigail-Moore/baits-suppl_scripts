
# Scripts for analyzing BUCKy data

R_masterscript_01.v5_REVISION.R:
This R script parses results from MrBayes analyses per locus.

R_masterscript_02.v4_REVISION.R:
This R script processes and visualizes results from MrBayes and BUCKy analyses.
It is designed to be executed after all the steps in R_masterscript_01.R

# Procedure for running the validation analyses using alignable transcriptome data:

1.  BLAST each of the transcriptomes against the same BLAST database used for the baits.

2.  Concatenate both the BLAST results and the transcriptome sequences into one file each.

3.  Parse the BLAST output, make separate files for each gene family, and add them to the backbone alignments for their gene families.
(Script 20170806_cary_trans_processing.py)

4.  Run tal_combiner.py (from the pipeline scripts) once to prune out all parts of the alignments that do not include baits sequences.
(This is because the transcriptome sequences are often longer than the sequences from the baits and it gets rid of the parts of the sequences
that cannot be reliably classified.)

5.  Run tal_combiner.py again to remove everything except the transcriptome sequences from the alignments.

6.  Take these alignments and slice them up randomly into 300 base long chunks (10 replicates were performed).
(Script 20170807_cary_trans_processing_pt2.py)

7.  Use fasta_renaming.py (from the pipeline scripts) to rename all of the sequences according to the format Individual_Name-Sequence_Name.
(They cannot be named in this format originally or tal_combiner.py will not keep the original sequence names.)

8.  Run them through Part II of the pipeline, using a revised version of tcontig_classif_master.py (tcontig_classif_master_validation.py) 
that calls a revised version of tbaits_intron_removal.py (tbaits_intron_removal_shortened.py) that takes the sequences as they are instead of
looking for output from BLASTing them against the exon database.

9.  Run them through Part III of the pipeline as normal.

10.  The output can be analyzed using the scripts ttest_results_processing.py (for the output of parts II and III) or ttest_results_processing_pt2only.py
(for the output of part II only).
