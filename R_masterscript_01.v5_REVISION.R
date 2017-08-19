
##################################################################################

### Jurriaan de Vos, 03 Aug 2017
### j.devos@kew.org

##################################################################################

## R_masterscript_01

##  This R script parses results from MrBayes analyses per locus and does 4 things

# 0. getting ready, loading libraries, etc.
# 1. it parses the raw MrBayes tree files:
#     - remove burnin and combine
#     - rename tips while scoring taxon sampling
#     - thin
# 2. it scores the posterior probability of monophyly 
#     - based on the thinned (200 trees) files
# 3. it prunes the trees to 1 random exemplar per family
#     - based on the 1500 postburnin trees, keeping 1 random exemplar
#     - keepers are renamed to their taxonomic family and saved (BUCKy input)
#     - keepers are further pruned by collapsing ACPT (more BUCKy input)
# 4. it sets up the required files for the BUCKy anaylses
#     - this requires doing mbsum and bucky as commands in a Terminal window.
#
# At the end of this script, and after running bucky, all results have been generated.
# Those will be visualized and processed using R_masterscript_02

#############################################################
##########   PART 0. getting ready
#############################################################

## required libraries
  library(ape)

## base working directory and helperfiles
  # change this to whatever works for you
  wd <- "D:/Dropbox/Paper_baits_2017/bucky/"

## read file that maps genera to their taxonomic family (following Nyffeler & Eggli 2010, Taxon)
  taxontableFileName <- paste(wd, "helperfiles/", "taxon_table_20170531.txt", sep="")

## make directories
  dir.create(paste(wd, "02_thinned_renamed_trees", sep=""))
  dir.create(paste(wd, "03_ready_pruned_trees", sep=""))
  dir.create(paste(wd, "04_mbsum_trees", sep=""))
  dir.create(paste(wd, "05_bucky_input", sep=""))
  dir.create(paste(wd, "06_bucky_output", sep=""))
  dir.create(paste(wd, "R_results", sep=""))

# MrBayes .t and .trprobs files are expected in the folder 01_raw_MrBayes_files
  # specifically, it expects two runs per locus
  # A file linking genera to family prefixes is expected in the folder helperfiles

# Thinning options
  nthinned <- 200 # thin pd to 200 trees
  burnin <- 0.25  # 25% burnin



#############################################################
##########   PART 1. parse trees, rename, thin
#############################################################


### read file names etc
  fileNames <- list.files(paste(wd, "01_raw_MrBayes_files/", sep=""))
  runs1 <- fileNames[grepl(".nex.run1.t", fileNames)]
  runs2 <- fileNames[grepl(".nex.run2.t", fileNames)]
  fileStems <- gsub(".nex.run1.t", "", runs1)

  
### parse each locus
  startTime <- Sys.time()
  for(locusnumber in 1:length(fileStems)) {
    #locusnumber <- 1  # debug
    cat(paste("\n", locusnumber, fileStems[locusnumber], "..."))
  
      
    ## read, prune burnin, combine runs, write files
    
    # read
    cat("\tRead 1...")
    trees1 <- readLines(paste(wd, "01_raw_MrBayes_files/", runs1[locusnumber], sep=""))
    cat("Read 2...")
    trees2 <- readLines(paste(wd, "01_raw_MrBayes_files/", runs2[locusnumber], sep=""))
    
    
    # check wether translation tables are identical
    cat("Check...")
    begin1 <- which(trees1 == "   translate")+1
    begin2 <- which(trees2 == "   translate")+1
    end1 <- grep("      [^ ]+ [^;]+;", trees1)  # regex for 5 spaces, a number, a word, a semicolon
    end2 <- grep("      [^ ]+ [^;]+;", trees2)
    if(FALSE %in% c(trees1[begin1:end1] == trees2[begin2:end2], length(begin1:end1)==length(begin2:end2))) {
      print("error; translation tables not identical - check files"); break
    } else {cat("ok! ")}
    
  
    # combine trees    
    cat("Combine...")
    # get lines that are trees post burnin
    trees_lines <- c(c(end1+1):c(grep("end;", trees1)-1))
    ntrees <- length(trees_lines)
    postBurninTreesIndex <- c(  c(round(ntrees*burnin)+1+trees_lines[1]):c(ntrees+trees_lines[1])  )[1:750]
    # combine
    trees <- c(trees1[1:end1], trees1[postBurninTreesIndex], trees2[postBurninTreesIndex], "end;")
    # fix numbering  
    trees <- gsub("tree gen.[^ ]+ =", "tree unnamed =", trees) 
  
      
    # write combined trees
    cat("Write...")
    fileToMake <- paste(wd, "02_thinned_renamed_trees/", fileStems[locusnumber], ".combined.trees", sep="")
    cat(trees, file=fileToMake, sep="\n")
  }
  endTime <- Sys.time()

### parse each combined-trees file to rename tips
  
  ## read taxon table that matches genera to family
  taxontable <- read.table(taxontableFileName, stringsAsFactors = FALSE, header=TRUE)
  # make result object to track locus properties
  fams <- c("OUT", "MOL", "MON", "BAS", "HAL", "DID", "TAL", "POR", "ANA", "CAC")
  locustable <- matrix(0, nrow=length(fileStems), ncol=c(length(fams)))
  colnames(locustable)  <- fams
  rownames(locustable)  <- fileStems

  
  ## read trees
  for(locusnumber in 1:length(fileStems)) {
    #locusnumber <- 1  # debug
    cat(paste("\n", locusnumber, fileStems[locusnumber], "...Reading..."))
    lines <- readLines(paste(wd, "02_thinned_renamed_trees/", fileStems[locusnumber], ".combined.trees", sep=""))
    
    ## rename genus-by-genus to add family prefix
    cat("Renaming...")
    for(i in 1:nrow(taxontable)) {
      #i <- 1  # debug
      # define replacement string
      was <- paste(taxontable$genus[i], "_", sep="")  # _ to make sure genera are whole words cf. Portulaca ria
      willbe <- paste(taxontable$family[i], was, sep="_")
      
      # count number of occurrences to track stats
      n <- sum(grepl(was, lines))
      thisRow <- fileStems[locusnumber]
      thisColumn <- which(colnames(locustable)==taxontable$family[i])
      locustable[thisRow, thisColumn]  <- locustable[thisRow, thisColumn] + n
      
      # replace
      lines <- gsub(was, willbe, lines, fixed = TRUE)    
    }
    cat("done.")    
    ## write, now that all replacements have been made
    cat(lines, file = paste(wd, "02_thinned_renamed_trees/", fileStems[locusnumber], ".renamed.trees", sep=""), sep="\n")
  }

  
  
###  thin tree files to specified number of trees (default 200)
  
  for(locusnumber in 1:length(fileStems)) {
    #locusnumber <- 1  # debug
    cat(paste("\n", locusnumber, fileStems[locusnumber], "...Reading..."))
    trees <- read.nexus(paste(wd, "02_thinned_renamed_trees/", fileStems[locusnumber], ".renamed.trees", sep=""))

    ## thin  renamed trees to target number to keep (nthinned)
    cat("Thinning...")
    tmp <- length(trees)  # trees number prior to thinning
    tmp <- seq(1, tmp, c(tmp/nthinned))  # equally spread trees will be taken 
    newThinnedTrees <- trees[tmp] # thin to target number

    ## write files
    cat("Writing...")
    fileToMake <- paste(wd, "02_thinned_renamed_trees/", fileStems[locusnumber], ".thinned.renamed.trees", sep="")
    write.nexus(newThinnedTrees, file=fileToMake)
    cat("Done.")
  }

# done

  ## Saving progress as checkpoint
  # tree files are already saved to files, but locustable and fileStems is required in memory for the next steps
  save.image(file=paste(wd, "R_results/checkpoint_1.Rdata", sep=""))





#################################################################
##########   PART 2. score monophyly of clades of interest
#################################################################


### parse each THINNED treeFile to check monophyly
  # monophyly is checked based on trying to root the tree by the tips from the clade to check
  # after removing outgroup sequences (i.e., everything not Portulacineae + Molluginaceae)
  
  # The clades that are checked are those that 

### define groups to check monophyly for 
  # names of clades of which we track monophyly: (defined below)
  # we do not check outgroup monophyly
  # the clades that are checked for monophyly are those suspected to be possibly supported
  tmp <- fams[fams!="OUT"]
  monogroups <- c(tmp, "AC", "PC", "TC", "AP", "ACP", "ACPT", "ACPTD", "ACPTH", "ACPTB", "ACPTM", "BH", "MH", "DH", "BM") # defined below
  monophylyCounters <- data.frame(matrix(0, nrow=length(fileStems), ncol = length(monogroups)))
  rownames(monophylyCounters) <- fileStems
  colnames(monophylyCounters) <- monogroups

  lociToExclude <- NULL  # object to track whether drop.tip encounters spurious loci

### run over loci
  for(locusnumber in 1:length(fileStems)) {
    #locusnumber <- 1  # debug
    loc <- fileStems[locusnumber]
    print(paste(locusnumber, loc))    
    
    ## read trees for locus and define summary objects
    fileToRead <- paste(wd, "02_thinned_renamed_trees/", loc, ".thinned.renamed.trees", sep="")
    trees <- read.nexus(file=fileToRead)
    ntrees <- length(trees)
    monophylyCounter <- rep(0, length(monogroups))  # to score posterior probability of groups being monophletic
    names(monophylyCounter) <- monogroups
    
    
    ## run over each tree
    for (i in 1:ntrees) {
      #i<-1  # debug
      
      ## read trees and drop outgroup
      tree <- trees[[i]]  # extract tree
      # drop outgroup (tips with "OUT_") from tree
      outtips <- tree$tip.label[grepl("OUT_", tree$tip.label)]
      if(length(outtips)+4 > length(tree$tip.label)) { # TRUE if tree after pruning has less than 4 taxa
                                                       # i.e. it doesnt contain relevant information
        if(i == 1){
          print(paste(locusnumber, fileStems[locusnumber], "has too few ingroup taxa; skipping..."))
          lociToExclude <- locusnumber
        }  
        next()
      }
      tree <- drop.tip(tree, outtips)
      
      
      ## make an object with subfamily of each tip
      # first 3 letters is always taxonomic family, 4th optional character for subfamily, _ otherwise
      tips <- substring(tree$tip.label, 1, 4)
      tips <- gsub("_", "", tips)
      
      
      ## track the monophyly of all the clades of interest (termed monogroups), except the outgroup
      for (monogroup in monogroups) {  
        # monogroup <- monogroups[1]
        
        # where are the tips of interest?
        # defining monogroups 
        if (monogroup %in% fams)     {searchphrase <- monogroup}
        if (monogroup == "AC")       {searchphrase <- c("CAC", "ANA")}
        if (monogroup == "PC")       {searchphrase <- c("CAC", "POR")}
        if (monogroup == "TC")       {searchphrase <- c("CAC", "TAL")}
        if (monogroup == "AP")       {searchphrase <- c("POR", "ANA")}
        if (monogroup == "ACP")      {searchphrase <- c("CAC", "POR", "ANA")}
        if (monogroup == "ACPT")     {searchphrase <- c("CAC", "POR", "ANA", "TAL")}
        if (monogroup == "ACPTD")    {searchphrase <- c("CAC", "POR", "ANA", "TAL", "DID")}
        if (monogroup == "ACPTH")    {searchphrase <- c("CAC", "POR", "ANA", "TAL", "HAL")}
        if (monogroup == "ACPTB")    {searchphrase <- c("CAC", "POR", "ANA", "TAL", "BAS")}
        if (monogroup == "ACPTM")    {searchphrase <- c("CAC", "POR", "ANA", "TAL", "MON")}
        if (monogroup == "ACPTDBHM") {searchphrase <- c("CAC", "POR", "ANA", "TAL", "DID", "BAS", "HAL", "MON")}
        if (monogroup == "BH")       {searchphrase <- c("BAS", "HAL")}
        if (monogroup == "MH")       {searchphrase <- c("MON", "HAL")}
        if (monogroup == "DH")       {searchphrase <- c("DID", "HAL")}
        if (monogroup == "BM")       {searchphrase <- c("BAS", "MON")}
        tipIndex <- which(tips %in% searchphrase)
        
        
        # process depending on how many tips were encountered
        if(length(tipIndex)>1) {
          # check if monophyletic by checking for an error when attempting to root by the monogroup
          # if so add to pp indicator
          # check monophyly on all tips (incl outgroup)
          tmptree <- try(root(tree, tipIndex), silent = TRUE)
          if("phylo" %in% class(tmptree)) {  # rooting succesfull; add contribution of tree to pp of clade
            monophylyCounter[monogroup] <- monophylyCounter[monogroup]+(1/ntrees)
          }
        }
        if(length(tipIndex)==1) {
          # score that its monophyletic  (even though it is not very informative...)
          monophylyCounter[monogroup] <- monophylyCounter[monogroup]+(1/ntrees)
        }
        if(length(tipIndex)<1) {
          # make a note in monophylycounter but don't do anything else
          monophylyCounter[monogroup] <- NA
        }
      }
    }
    ## save monophylyCounter object for this locus
    monophylyCounters[locusnumber, ] <- monophylyCounter
  } 
  # done with all loci    

  
### save image as checkpoint

  save.image(file=paste(wd, "R_results/checkpoint_2.Rdata", sep=""))
  



#################################################################
##########   PART 3. prune families and ACPT to 1 random exemplar
#################################################################


# pruning to a random exemplar per tree; given a lot of trees used for bucky per locus, 
# this accounts for uncertainty in whether a family is monophyletic or not

  collapse_to_one_tip_per_family <- function(tree) {
    # helperfunction to prune a tree to one *random* exemplar per family, which are 
    # taken to be indicated by the first three letters of the tip.label
    # execution with help of subhelperfunction because sample() is funky:
    #  if there's only one taxon, sample(which(..), 1) works counterintuitively: 
    #  it doesnt pick the only one element in which(..), 
    #  but it picks a number from the range 1:which(..)  # (thanks, R)
    resample <- function(x, ...) {
      if(length(x > 0)) {
        x[sample.int(length(x), ...)]
      } else {x}
    }
    
    # the families present
    fams_among_tips <- unique(substr(tree$tip.label, 1, 3))
    # an object to store selected tips
    selected_tips <- NULL
    
    for(fam in fams_among_tips) {
      tips <- grep(fam, tree$tip.label, value=TRUE)
      selected_tips <- c(selected_tips, resample(tips,1))
    }
    
    # prune to keep selected_tips
    tree2 <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% selected_tips])
    
    # rename all tips
    tree2$tip.label <- substr(tree2$tip.label, 1, 3)
    
    return(tree2)
  }  



  collapse_ACPT_to_one_tip <- function(tree) {
    # helperfunction to prune a tree to one *random* exemplar for ACPT, which is taken to be
    # indicated by CAC ANA POR and TAL
    # execution with help of subhelperfunction because sample() is funky. (see above)
    
    # input is supposed to be a renamed, pruned tree from function collapse_to_one_tip_per_family()
    
    # helperfunction
    resample <- function(x, ...) {
      if(length(x > 0)) {
        x[sample.int(length(x), ...)]
      } else {x}
    }
    
    ACPTindex <- which(tree$tip.label %in% c("CAC", "ANA", "POR", "TAL"))
    notACPTindex <- which(!tree$tip.label %in% c("CAC", "ANA", "POR", "TAL"))

    # exception for if there's no non-ACPT tips
    # then we cannot drop; gotta skip
    if(length(notACPTindex) < 1) {
      tree2 <- NULL  
    } else { # there will be a tree returned
      
      # drop if needed
      if(length(ACPTindex) > 1) { # we need to drop something
        to_keep <- resample(ACPTindex, 1)
        to_drop <- ACPTindex[!ACPTindex %in% to_keep]
        tree2 <- drop.tip(tree, tree$tip.label[to_drop])
      } else {                    # we dont need to drop something
        tree2 <- tree
      }
      # rename remaining ACPT tip
      # rename
      tree2$tip.label[tree2$tip.label %in% c("CAC", "ANA", "POR", "TAL")] <- "ACPT"
    }
    return(tree2)
  }  

  
  ### use the helper function to prune to one random tip per family

  for(locusnumber in 1:length(fileStems)) {
    #locusnumber <- 1
    # load trees; collapse; rename; save
    loc <- fileStems[locusnumber]
    cat(paste("\n", locusnumber, loc, "\tReading..."))
    
    # load
    trees <- read.nexus(paste(wd, "02_thinned_renamed_trees/", fileStems[locusnumber], ".renamed.trees", sep=""))
    
    # prune
    cat("Pruning...")
    trees2 <- lapply(trees, collapse_to_one_tip_per_family)

    # save
    cat("Saving...")
    fileToMake <- paste(wd, "03_ready_pruned_trees/", fileStems[locusnumber], ".pruned.trees", sep="")
    write.nexus(trees2, file=fileToMake)  ## this is the slow step
    cat("done!")
  }


  ### use the collapse ACPT function to prune to one random tip for ACPT
  
  for(locusnumber in 1:length(fileStems)) {
    #locusnumber <- 1
    # load trees; collapse; rename; save
    loc <- fileStems[locusnumber]
    cat(paste("\n", locusnumber, loc, "\tReading..."))
    
    # load
    trees <- read.nexus(paste(wd, "03_ready_pruned_trees/", fileStems[locusnumber], ".pruned.trees", sep=""))
    
    # prune
    cat("Pruning...")
    trees2 <- lapply(trees, collapse_ACPT_to_one_tip)  # returns trees if pruning made sense, FALSE otherwise
    
    # check if successfull
    if(inherits(trees2[[1]], "NULL")) {
      cat("too few non-ACPT; skipping...")
    } else {
      # save
      cat("Saving...")
      fileToMake <- paste(wd, "03_ready_pruned_trees/", fileStems[locusnumber], ".pruned_noACPT.trees", sep="")
      write.nexus(trees2, file=fileToMake)
    }
    cat("done!")
  }


  ##### Wrapping up, two more things to do:  

  ## 1. join monophylyCounter object and taxon sampling object and save to txt table
  
  tmp1 <- locustable
  colnames(tmp1) <- paste("n", colnames(locustable), sep="_")
  tmp2 <- monophylyCounters
  pps <- sapply(tmp2, is.numeric) & !sapply(tmp2, is.integer)  # round to two digits to avoid float errors
  tmp2[, pps] <- sapply(tmp2[, pps], round, digits=2)
  colnames(tmp2) <- paste("pp", colnames(monophylyCounters), sep="_")
  
  allLocusResults <- cbind(tmp1, tmp2[rownames(tmp1), ])

    # workaround to make header of file align to columns (since write.table(..., row.names=T))
  tmp <- allLocusResults
  colnames(tmp)[1] <- paste("locus\t", colnames(tmp)[1], sep="")
  write.table(tmp, file=paste(wd, "R_results/allLocusResults.txt", sep=""), quote=FALSE, sep="\t")
  
  
  ## 2. rename all the trees in the .../03.../ folder so that bucky understands them.

  files <- list.files(paste(wd, "03_ready_pruned_trees", sep=""))
  for(file in files) {
    # thanks https://gist.github.com/mages/1544009
    print(file)
    fileToEdit <- paste(wd, "03_ready_pruned_trees/", file, sep="")
    x <- readLines(fileToEdit)
    y <- gsub("TREE * ", "tree ", x, fixed = TRUE)    
    cat(y, file = fileToEdit, sep="\n")
  }
  
  
  # done; save checkpoint
  save.image(file=paste(wd, "R_results/checkpoint_3.Rdata", sep=""))
  
  
  #################################################################
  ##########   PART 4. prepare, set up and run BUCKy analyses
  #################################################################
  
  
  #### step 1: run mbsum to parse MrBayes-like output trees to bucky's format from the .../03_ready_pruned_trees/ folder
  #### This bit requires running BUCKy from the Terminal: 
  # for FILE in *pruned.trees; 
  #   do mbsum -n 0 -o ../04_mbsum_trees/$FILE.mbsum $FILE; 
  # done
  
  
  #### step 2:, define 3 BUCKy analyses using a helper function:
  # - one with all ACPT families
  # - one with all Portullugo families
  # - one with ACPT collapsed and all other Portollugo families
  # a helper function copies the selected loci to a new folder
    

  copyBuckyInput <- function(targetLoci, dirName, suffix) {
    # helperfunction to write locus mbsum files to the location of input loci
    # targetLoci are the names (as in fileStems and rownames(allLocusResults))
    # dirname is the name of the subdirectory where the files are to be copied to
    # suffix is the identifier suffix of the locus name (e.g. ".pruned.trees.mbsum")

    # make target directory
    targetDir <- paste(wd, "05_bucky_input/", dirName, "/", sep="")
    dir.create(targetDir)
    # copy desired files
    source_files <- paste(wd, "04_mbsum_trees/", targetLoci, suffix, sep="")
    for(i in 1:length(source_files)) {
      file.copy(source_files[i], targetDir)
    }
  }  

  
  ####  write bucky input files using helperfunction
  x <- allLocusResults  # unless you're not into the whole brevity thing

  ## Portullugo analysis
  targetLoci <- rownames(x)[
      x$n_MOL > 0 &
      x$n_MON > 0 &
      x$n_BAS > 0 &
      x$n_HAL > 0 &
      x$n_DID > 0 &
      x$n_TAL > 0 &
      x$n_POR > 0 &
      x$n_ANA > 0 &
      x$n_CAC > 0
    ]
  copyBuckyInput(targetLoci, "allfams", ".pruned.trees.mbsum")
 
  ## non-ACPT analysis
  targetLoci <- rownames(x)[
      x$n_MON > 0 &
      x$n_BAS > 0 &
      x$n_HAL > 0 &
      x$n_DID > 0 &
     (x$n_TAL + x$n_POR + x$n_ANA + x$n_CAC) > 0 &
      x$n_MOL > 0
    ]
  copyBuckyInput(targetLoci, "nonACPT", ".pruned_noACPT.trees.mbsum")

  ## ACPT analysis
  targetLoci <- rownames(x)[
    x$n_TAL > 0 &
    x$n_POR > 0 &
    x$n_ANA > 0 &
    x$n_CAC > 0
    ]
  copyBuckyInput(targetLoci, "ACPT", ".pruned.trees.mbsum")
  

  #### step 3: run BUCKy from the terminal
  # run it from the folder 05_bucky_input:
  # for DIR in *; 
  #   do bucky --calculate-pairs -k 4 -o ../06_bucky_output/$DIR $DIR/*;
  # done
  
  # this does a default analysis but with 4 runs and calculates the probability 
  # for each pair of loci that those two loci suppotr the same tree

  
### all done.
  
### now visualize and analyze results using R_masterscript_02
  
  

  
  
  
  
  
  
  
  
  
  
  
  