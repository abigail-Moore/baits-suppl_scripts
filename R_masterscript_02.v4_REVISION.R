
##################################################################################

### Jurriaan de Vos, 03 AUG 2017
### j.devos@kew.org

##################################################################################

## R_masterscript_02

##  This R script processes and visualizes results from MrBayes and BUCKy analyses.
##  It is designed to be executed after all the steps in R_masterscript_01.R

##  It does several things:

# 0. getting ready, loading libraries, etc.
# 1. plotting MrBayes MAP trees with nice intuitive colours
# 2. plotting histograms of the support values across gene trees for specific nodes
# 3. plotting bucky results from the calculate-pairs option in the form of custom heatmaps


#############################################################
##########   PART 0. getting ready
#############################################################

## required libraries
  library(ape)
  library(RColorBrewer)

## base working directory and helperfiles
# change this to whatever you used in R_masterscript_01
  wd <- "D:/Dropbox/Paper_baits_2017/bucky/"

# get the necessary objects
  load(paste(wd, "R_results/checkpoint_3.Rdata", sep=""))
  
  
  
    
################################################################
##########   PART 1. plotting MrBayes Maximum A Posteriori trees
################################################################
  
  # ##  the goal is to plot every gene tree (in the form of a mrbayes Maximum a-posteriori probability tree, MAP)
  #  #  with attractive colours that aid interpretation, and reporting posterior probabilities of relevant nodes
  # 
  # ## helper function that does all the hard work
  # 
  # plotRootedTreeWithColors <- function(tree, label="", pps, q=7, allPresentCriterion=FALSE, ...) {
  #   ## this function plots a tree with nice colors
  #   ## first it tries to root it in various ways
  #   ## and then plots it with nice colors based on subfamily affinity
  #   ## arguments:
  #    # tree: a tree (e.g., MAP), 
  #    # label: a label that is plotted in the lower left corner
  #    # pps: a line of the allLocusResults files with pp of clades to be plotted
  #    # q: a plotting parameter q, where 1/q sets the fraction of space to the left of pp numbers
  #    # allPresentCriterion: a switch, which determines whether pp of clades to be reported require
  #    #   at least one representative of each family
  #    # ... for additional input to text()
  # 
  #   ## Step 1, root the tree using the following criteria in order of preference:
  #   # 1. the mrca of the present model taxa 
  #   #    oryza, vitis, solanum, arabidopsis, glycine
  #   # 2. the mrca of all taxa with prefix OUT or MOL *if monophyletic*
  #   # 3. a random taxon with prefix OUT
  #   # 4. a random taxon with prefix MOL
  #   # 5. a random taxon with prefix MON
  #   # 6. a random taxon with prefix BAS
  #   # 7. a random taxon with prefix DID
  #   # 8. a random taxon with prefix TAL
  #   # 9. a random taxon with prefix POR
  #   #10. a random taxon with prefix ANA
  #   #11. a random taxon with prefix CAC
  #   
  #   rootSuccess <- 0 # tracker
  #   
  #   # criterion 1
  #   models <- c("Oryza", "Vitis", "Solanum", "Arabidopsis", "Glycine")
  #   focalTips <- sapply(models, function(x) {tree$tip.label[grepl(x, tree$tip.label)]})
  #   focalTips <- unlist(focalTips[!is.na(focalTips)])
  #   if(length(focalTips) > 0) { # root criterion 1 is met
  #     if(length(focalTips)==1) {
  #       tree_root <- root(tree, focalTips)
  #     } else {tree_root <- root(tree, node=getMRCA(tree, focalTips))}
  #     rootSuccess <- 1
  #   }
  #   
  #   # criterion 2
  #   focalTips <- sapply(c("OUT", "MOL"), function(x) {tree$tip.label[grepl(x, tree$tip.label)]})
  #   focalTips <- unlist(focalTips)
  #   tmp <- try(root(tree, focalTips), silent = TRUE)
  #   if(rootSuccess==0 & "phylo" %in% is(tmp)) {
  #     if(length(focalTips)==1) {
  #       tree_root <- root(tree, focalTips)
  #     } else {tree_root <- root(tree, node=getMRCA(tree, focalTips))}
  #     rootSuccess <- 1
  #   }  # otherwise, root threw an error that indicated outgroup wasnt monophyletic
  #   
  #   # criteria 3-11
  #   roots <- c("OUT", "MOL", "MON", "BAS", "DID", "TAL", "POR", "ANA", "CAC")
  #   for(i in roots) {
  #     if(rootSuccess==0 & (TRUE %in% grepl(i, tree$tip.label))) {
  #       tree_root <- root(tree, sample(grep(i, tree$tip.label, value=TRUE), 1))
  #       rootSuccess <- 1
  #     }
  #   }
  #   
  #   ## ladderize tree to make nice and aid interpretation
  #   tree_root <- ladderize(tree_root)
  #   
  #   ## make nice colors to plot
  #   # there's 13 possible subfamilies; grab 9 cool colors and a 10th; replace gray and yellow for enhanced aesthetics
  #   cols <- c(brewer.pal(9, "Set1"), "black")
  #   cols[9] <- brewer.pal(9, "Greens")[7]
  #   cols[3] <- brewer.pal(9, "Greens")[4]
  #   cols[6] <- brewer.pal(9, "Blues")[4]
  #   
  #   tipcolIndex <- tree_root$tip.label  
  #   tipcolIndex[grepl("OUT", tipcolIndex)] <- 10
  #   tipcolIndex[grepl("MOL", tipcolIndex)] <- 1
  #   tipcolIndex[grepl("MON", tipcolIndex)] <- 2
  #   tipcolIndex[grepl("BAS", tipcolIndex)] <- 3
  #   tipcolIndex[grepl("HAL", tipcolIndex)] <- 4
  #   tipcolIndex[grepl("DID", tipcolIndex)] <- 5
  #   tipcolIndex[grepl("TAL", tipcolIndex)] <- 6
  #   tipcolIndex[grepl("POR", tipcolIndex)] <- 7
  #   tipcolIndex[grepl("ANA", tipcolIndex)] <- 8
  #   tipcolIndex[grepl("CAC", tipcolIndex)] <- 9
  #   tipcolIndex <- as.numeric(tipcolIndex)
  #   tipcols <- cols[tipcolIndex]  
  #   
  #   ## check if any names need updating
  #   tree$tip.label[tree$tip.label == "MOL_Mollugo_WA_endemic_137"]  <- "MOL_Trigastrotheca_molluginea_137"
  #   tree$tip.label[tree$tip.label == "MOL_Mollugo_pentafilis_148"]  <- "MOL_Trigastrotheca_stricta_148"
  #   tree$tip.label[tree$tip.label == "CAC_Pereskia_saccharosa_67"]  <- "CAC_Pereskia_sacharosa_67"
  #   tree$tip.label[tree$tip.label == "POR_Portulaca_sp_209"]        <- "POR_Portulaca_cf_perennis_209"
  #   tree$tip.label[tree$tip.label == "POR_Portulaca_phillippi_213"] <- "POR_Portulaca_cf_perennis_213"
  #   tree$tip.label[tree$tip.label == "POR_Portulaca_sp_clay_pan_187"] <- "POR_Portulaca_cf_filifolia_187"
  #   tree$tip.label[tree$tip.label == "POR_Portulaca_sp_184"]        <- "POR_Portulaca_cf_digyna_184"
  #   tree$tip.label[tree$tip.label == "POR_Sedopsis_filsonii_186"]   <- "POR_Portulaca_filsonii_186"
  #   tree$tip.label[tree$tip.label == "CAC_Blossfeldia_204"]         <- "CAC_Blossfeldia_liliputana_204"
  #   tree$tip.label[tree$tip.label == "CAC_Maihuenia_poepiggii_205"] <- "CAC_Maihuenia_poeppigii_205"
  #   tree$tip.label[tree$tip.label == "TAL_Talinum_arnotii_4"]       <- "TAL_Talinum_arnottii_4"
  # 
  #   ## plot the tree
  #   plot(tree_root, no.margin=TRUE, cex=0.6, align.tip.label=T, tip.color=tipcols)
  #   
  #   ## plot a selection of posterior probabilies
  #   # what clades to we report? These are the monogroups defined in masterscript 1 except the exotic TC
  #   what <- rev(c("pp_AP", "pp_AC", "pp_PC", "pp_ACP", "pp_ACPT", "pp_ACPTD", "pp_ACPTH", "pp_ACPTB", "pp_ACPTM", "pp_BH", "pp_MH", "pp_DH", "pp_BM"))
  # 
  #   # how do we report?
  #   if(allPresentCriterion==TRUE) {
  #     # then, we report a NA when not at least one member of each family is present
  #     if( !(pps$n_ANA > 0 & pps$n_POR > 0)) {pps$pp_AP <- NA}
  #     if( !(pps$n_ANA > 0 & pps$n_CAC > 0)) {pps$pp_AC <- NA}
  #     if( !(pps$n_POR > 0 & pps$n_CAC > 0)) {pps$pp_PC <- NA}
  #     if( !(pps$n_ANA > 0 & pps$n_CAC > 0 & pps$n_POR > 0))  {pps$pp_ACP <- NA}
  #     if( !(pps$n_ANA > 0 & pps$n_CAC > 0 & pps$n_POR > 0 & pps$n_TAL > 0))  {pps$pp_ACPT <- NA}
  #     
  #     if( is.na(pps$pp_ACPT)| !(pps$n_DID > 0)) {pps$pp_ACPTD <- NA}
  #     if( is.na(pps$pp_ACPT)| !(pps$n_HAL > 0)) {pps$pp_ACPTH <- NA}
  #     if( is.na(pps$pp_ACPT)| !(pps$n_BAS > 0)) {pps$pp_ACPTB <- NA}
  #     if( is.na(pps$pp_ACPT)| !(pps$n_MON > 0)) {pps$pp_ACPTM <- NA}
  #     
  #     if( !(pps$n_BAS > 0 & pps$n_HAL > 0)) {pps$pp_BH <- NA}
  #     if( !(pps$n_MON > 0 & pps$n_HAL > 0)) {pps$pp_MH <- NA}
  #     if( !(pps$n_DID > 0 & pps$n_HAL > 0)) {pps$pp_DM <- NA}
  #     if( !(pps$n_BAS > 0 & pps$n_MON > 0)) {pps$pp_BM <- NA}
  #     
  #   }
  #       
  #   # the y axis depends on the number of tips: use the lower third for the pps
  #   y <- length(tree_root$tip.label)/3.3  # highest pp line
  #   yseq <- seq(1, y, length.out=c(length(what)+1))
  #   
  #   # for x axis, extract current coords of plot
  #   lims <- par("usr")
  #   x <- ((lims[2]-lims[1]) / q) + lims[1] # # plot to 1/qth of the range to the right
  #   
  #   for(i in 1:length(what)) {
  #     lab <- what[i]
  #     text(x=0, y=yseq[i+1], labels=lab, pos=4, cex=0.8) # +1 to skip 0 when plotting
  #     lab <- sprintf("%.2f", pps[what[i]])
  #     text(x=x, y=yseq[i+1], labels=lab, pos=4, cex=0.8) # +1 to skip 0 when plotting
  #   }  
  #   # plot the label (e.g. a locus name)
  #   text(x=0, y=1, labels=label, pos=4, cex=0.8, font=2)
  #   
  #   ## done plotting
  # }
  # 
  # ## use the plotting helperfunction to plot MAP trees directly read from MrBayes file  
  # # the required MAP files are produced by MrBayes and expected in the
  # # .../01_raw_MrBayes_files/ folder
  # 
  # tmp <- list.files(paste(wd, "01_raw_MrBayes_files", sep=""))
  # treeFiles <- tmp[grepl("trprobs", tmp)]
  # 
  # pdf(paste(wd, "R_results/locusMAPTreePlots.pdf", sep=""), width=8.5, height=11)
  # par(mfrow=c(3, 2))
  # for(i in 1:length(treeFiles)) {
  #   #for(i in 1:6) {
  #   #i<-127
  #   #i<-3
  #   pd <- treeFiles[i]
  #   print(paste(i, pd))
  #   
  #   # read MAP tree
  #   tree <- read.nexus(paste(wd, "01_raw_MrBayes_files/", pd, sep=""))
  #   
  #   if(inherits(tree, "multiPhylo")) {tree <- tree[[1]]}
  #   if(inherits(tree, "phylo")) {tree <- tree}
  #   
  #   # rename tips
  #   tips <- tree$tip.label
  #   for(j in 1:length(tips)) {
  #     # extract genus name
  #     genus <- strsplit(tips[j], "_")[[1]][1]
  #     # match to family
  #     family <- taxontable[which(taxontable[,"genus"]==genus), 2]
  #     # rename
  #     tips[j] <- paste(family, tips[j], sep="_")
  #   }
  #   tree$tip.label <- tips  # reinsert
  # 
  #   # make label    
  #   loc <- rownames(locustable)[i]
  #   loc <- sub("c2p1pgts2_", "", loc)
  #   loc <- sub("_mb", "", loc)
  #   label <- paste("Locus ", i, ": ", loc, sep="")
  #   
  #   # provide posterior probabilies
  #   pps <- allLocusResults[rownames(locustable)[i], ]
  #   
  #   # plot
  #   plotRootedTreeWithColors(tree, label, pps, allPresentCriterion=TRUE, q=5.5)
  # }
  # dev.off()  
  

  ################################################################
  ##########   PART 2. plotting support histograms per locus
  ################################################################
  
  # this plot represents the posterior probability support for competing sistergroup 
  # relations among the focal loci
  # Panels indicate support for relations among ACP
  
  # There are two versions of the plots, one for the positions within ACPT, and one
  #   for the position of Halophytaceae within Portulacineae
  
  # helperfunction to draw stats on each graph  
  getlocsupport <- function(x) {    # function to extract summary of support per locus
    # from a numeric vector
    x1<- x[!is.na(x)]
    nosup <- sprintf("%.1f", length(x1[x1<=0.05])/length(x1)*100)
    sup   <- sprintf("%.1f", length(x1[x1>=0.95])/length(x1)*100)
    result <- c(paste("pp <= 0.05:  ", nosup, "%", sep=""),  # % of loci without support
                paste("pp >= 0.95:  ", sup, "%", sep=""))    # % of loci with high support
    return(result)
  }

  # helperfunction to plot a histogram
  histplotter <- function(what, title, drawSupport=TRUE, ...) {
    hist(what, breaks=20, main="", xlim=c(0,1), xlab="", ylab="", xaxt="n", ...)
    axis(side=1, at=c(seq(0,1,0.2)), labels=c("0.0", "", "", "", "", "1.0"))
    text(x=0.05, y=c(145*ymax/150), title, pos=4, cex=cexText)
    mtext(side=1, line=2, "Posterior probability\nof focal node", cex=0.6)
    if(drawSupport) {
      x <- getlocsupport(what)
      text(x=0.05, y=c(135*ymax/150), x[1], pos=4, cex=cexText)
      text(x=0.05, y=c(125*ymax/150), x[2], pos=4, cex=cexText)
    }
  }  
  
  ###  version for within-ACPT relations
  pps <- allLocusResults

  pdf(file=paste(wd, "R_results/histo_supportedTrees_ACPT_simple_REVISION.pdf", sep=""), width=7.5*(3/4), height=2)
  par(mfrow=c(1,3), mar=c(4, 2.8, 1.5, 0.4))
  
  criterion1 <- rep(TRUE, nrow(pps))  # all TRUE; plotting all loci
  criterion2 <- pps$pp_ACPT >= 0.05   # to play around to exclude loci that do not support the ACPT clade, which is an uncontroversial clade
  criterion2[is.na(criterion2)] <- FALSE # only TRUE when ACPT occasionally present in posterior distro
  # sum(criterion)                               #  70 loci with support for ACPT
  # sum(pps$pp_ACPT[!is.na(pps$pp_ACPT)] < 0.05) #  73 loci with no ACPT support
  # nrow(pps)                                    # 200 loci total

  cexText <- 0.7766
  ymax <- 125
  # 

  histplotter(pps$pp_AP[criterion1], "(POR,ANA)", ylim=c(0,ymax), yaxt="n")
  axis(2, at=c(0,25,50,75,100, 125), tcl=-0.5, labels=c("0", "", "50", "", "100", ""))

  histplotter(pps$pp_AC[criterion1], "(CAC,ANA)", ylim=c(0,ymax), yaxt="n")
  axis(2, at=c(0,25,50,75,100, 125), tcl=-0.5, labels=c("0", "", "50", "", "100", ""))

  histplotter(pps$pp_PC[criterion1], "(CAC,POR)", ylim=c(0,ymax), yaxt="n")
  axis(2, at=c(0,25,50,75,100, 125), tcl=-0.5, labels=c("0", "", "50", "", "100", ""))

  dev.off()
  

  # ###  version for Halophytaceae-specific focus
  # pdf(file=paste(wd, "R_results/histo_supportedTrees_nonACPT.pdf", sep=""), width=7.5/4*3, height=2)
  # par(mfrow=c(1,3), mar=c(4, 2.8, 1.5, 0.4))
  # 
  # pps <- allLocusResults  
  # 
  # cexText <- 0.7766
  # ymax <- 150
  # 
  # histplotter(pps$pp_BH, "(HAL, BAS)", ylim=c(0,ymax), yaxt="n")
  # axis(2, at=c(0,25,50,75,100, 125, 150), tcl=-0.5, labels=c("0", "", "50", "", "100", "", "150"))
  # 
  # histplotter(pps$pp_MH, "(HAL,MON)", ylim=c(0,ymax), yaxt="n")
  # axis(2, at=c(0,25,50,75,100, 125, 150), tcl=-0.5, labels=c("0", "", "50", "", "100", "", "150"))
  # 
  # histplotter(pps$pp_ACPTH, "(HAL,ACPT)", ylim=c(0,ymax), yaxt="n")
  # axis(2, at=c(0,25,50,75,100, 125, 150), tcl=-0.5, labels=c("0", "", "50", "", "100", "", "150"))
  # 
  # dev.off()
  
  
################################################################
##########   PART 3. plotting BUCKy results
################################################################
  
  # aim: a heatmap figure with side panels that contain the competing topologies
  # to this end we need two files 
  # - a .concordance file with the CF support values
  # - a .pairs file with the pp that loci have the same tree
  
  # outline:  
  # Section 1. define helperfunction to extract relevant information from a .concordance file
  # Section 2. using the function to get the CF support values
  # Section 3. draw the heatmaps
  
  # Then, three heatmaps are plotted. Basically the same chunck of code with
  #   minor adjustments to the layout() and monophyly statistics reported

  
  ############  SECTION 1  ###########    
  ### define helper function to read bucky result file
  
  parse_buckyResultFile <- function(buckyResultFile, root=NULL) {
    ## parses a .concordance file from a bucky analyses and extracts results
    #  including Concordance Factors (CFs)
    #  returns: 
    #  - the primary concordance tree with sample-wide CF mean in a format that is 
    #    plotable with plot(p); drawSupportOnEdges(p$node.label)
    #  - the CFs for all splits with pp>0.05
    #  - the translation table to link tips to tip names (for cross-referencing / checking)
    
    #buckyResultFile <- buckyResultFiles[1]  # debug
    
    ## read tree by parsing the string at the line one below the header    
    x <- readLines(paste(wd, "06_bucky_output/", buckyResultFile, sep=""))
    
    ## get primary concorance tree
    headerLine <- grep("Primary Concordance Tree with Sample Concordance Factors:", x)
    pc_tree_raw <- x[headerLine+1]
    p <- read.tree(text=pc_tree_raw)
    
    ## extract tip translation table
    nTips <- length(p$tip.label)
    translationTable <- matrix(NA, nrow=nTips, ncol=2)
    for(i in 1:nrow(translationTable)) {
      translationTable[i, 1] <- gsub(" ", "", substr(x[i+1], 1, 3))
      translationTable[i, 2] <- gsub("[,;]", "", substr(x[i+1], 4, 99))  #99 for simply to the end
    }  
    translationTable <- data.frame(translationTable, stringsAsFactors=F)    
    
    ## get all splits (in the primary tree or not)
    headerLine1 <- grep("Splits in the Primary Concordance Tree", x)
    headerLine2 <- grep("Splits NOT in the Primary Concordance Tree but with estimated CF > 0.050", x)
    endLine <- grep("Average SD of mean sample-wide CF", x)
    splits_raw <- x[ c((headerLine1+1):(headerLine2-2), (headerLine2+1):(endLine-2))  ]
    splits_raw <- gsub("\t.*", "", splits_raw)    
    splits_raw <- strsplit(splits_raw, " ")
    nsplits <- length(splits_raw)
    splits <- data.frame(splitstring=rep(NA, nsplits), CF_sample=NA, CF_genome=NA, FewestMembers=NA)
    for(i in 1:nsplits) {splits[i, 1:3] <- splits_raw[[i]][1:3]}
    
    ## get the members on either side of the split, translate, sort, and merge
    tmp <- gsub("[{}]", "", splits$splitstring)
    for(i in 1:length(tmp)) {
      t <- tmp[i]
      t1 <- strsplit(t, "|", fixed=T)[[1]]
      t2 <- sapply(t1, function(x) gsub("[^,]", "", x))
      t3 <- sapply(t2, nchar)
      FewestMembers <- strsplit(t1[which.min(t3)], ",")[[1]]
      MostMembers <- strsplit(t1[which.max(t3)], ",")[[1]]
      
      # translate
      for(oldTip in translationTable[, 1]) {
        FewestMembers[FewestMembers %in% oldTip] <- translationTable[oldTip, 2]
        MostMembers[MostMembers %in% oldTip] <- translationTable[oldTip, 2]
      }
      # sort and merge
      splits[i, "FewestMembers"] <- paste(sort(FewestMembers), collapse=",")
      splits[i, "MostMembers"] <- paste(sort(MostMembers), collapse=",")
      
    }
    
    ## employ translation table on tree
    for(oldTip in translationTable[, 1]) {
      p$tip.label[p$tip.label==oldTip] <- translationTable[oldTip, 2]
    }
    
    ##  extract the support factor as node label
    #  the node information is stored as branch length.  this is now transfered to a $node.label:
    #  tree$node.label[i] is the length of the branch ending in i
    # so that is the branch that has as property that the second element of $edge is i
    m <- p$Nnode
    n <- length(p$tip.label)
    p$node.label<- rep(NA, m)  # There is no root node so m+2=n
    for(i in (n+1):(n+m))  { # for the numbers of all internal nodes
      r <- which(p$edge[,2] == i)  # node i results from this edge
      l <- p$edge.length[r]        # which has this length
      t <- i-(n)                   # which should go to this node label (nodenr n+2  is labelnr 2)
      if(length(l) > 0) p$node.label[t] <- l   # assign if the node comes with a label
    }
    
    ## root if a root tip is specified in the function call
    if(inherits(root, "character")) {
      p <- compute.brlen(root(p, root, resolve.root=T, edgelabel=T))
    }
    
    ## return the tree, all the splits, and the translation table
    return(list(phylo=p, splits=splits, translation=translationTable))
  }

  ############  SECTION 2  ###########    
  ### use helperfunctions to get CF info and plot trees

  # pdf(paste(wd, "R_results/BUCKy_PrimaryConcordancetrees.pdf", sep=""), width=3, height=4.25)
  # 
  # par(mfrow=c(3,1), mar=c(2,2,2,2))
  # 
  # res1 <- parse_buckyResultFile("nonACPT.concordance", root="MOL")
  # plot(p <- res1[[1]], label.offset= 0.2)    
  # drawSupportOnEdges(p$node.label, frame="none", adj=c(0.5,-0.5), cex=0.85)
  # 
  # res2 <- parse_buckyResultFile("ACPT.concordance", root = "TAL")
  # plot(p <- res2[[1]], label.offset= 0.2)    
  # drawSupportOnEdges(p$node.label, frame="none", adj=c(0.5,-0.5), cex=0.85)
  # 
  # res3 <- parse_buckyResultFile("allfams.concordance", root = "MOL")
  # plot(p <- res3[[1]], label.offset= 0.2)    
  # drawSupportOnEdges(p$node.label, frame="none", adj=c(0.5,-0.5), cex=0.85)
  # 
  # dev.off()
  
  
  # ###  report all the splits that BUCKy found into a table
  # 
  # splits <- rbind(cbind(analysis="allfams", res3$splits), 
  #       cbind(analysis="nonACPT", res1$splits),
  #       cbind(analysis="ACPT", res2$splits))
  # write.table(splits, file=paste(wd, "R_results/BuckySplitsTable.txt", sep=""), 
  #             sep="\t", col.names=T, row.names=F, quote=F)
  

  ############  SECTION 3  ###########    
  ### plotting the heatmap
  
  # name of the analysis we're plotting
  analysis <- "ACPT"  
  pdf(paste(wd, "R_results/BUCKY_heatmap_", analysis, "REVISION.pdf", sep=""), width=8.5, height=8.5)
  
  ######## 1. parse the input from files into R object, using the .pairs and .input files

  ### get the "pairs" files for the ACPT analysis
  # read input
  raw <- readLines(paste(wd,"06_bucky_output/", analysis, ".pairs", sep=""))
  # make single spaces as separators
  raw <- gsub(" +", " ", raw)
  # split; each row with locus information becomes a list element
  raw <- strsplit(raw, " ")
  # now its a list; make it a matrix
  # input contained 143 loci
  nlocus <- length(raw)  ## right??
  mat <- matrix(NA, nrow=nlocus, ncol=nlocus)
  for(i in 1:nlocus){
    # fill in matrix by row; for each list element start at the third element
    # because the first is empty and the second one is the locus number
    mat[i, ] <- as.numeric(raw[[i]][seq(3, nlocus+2, 1)])
  }
  
  ### make object linking locusnumbers in the matrix to locusnames
  # note that bucky counts from 0 and R counts from 1
  # read names from input file
  raw <- readLines(paste(wd,"06_bucky_output/", analysis, ".input", sep=""))
  raw<- raw[3:c(2+nlocus)]
  # extract locus names
  raw <- sub(paste(" *[0-9]+ ", analysis, "/", sep=""), "", raw)
  pairsLocusNames <- sub(".pruned.trees.mbsum", "", raw)
  
  ### use locus names for rows and columns of matrix
  rownames(mat) <- pairsLocusNames
  colnames(mat) <- pairsLocusNames
  
  #### 2. draw heatmap (relevant bits extracted from heatmap() )
  
  ## Define the plot layout
  tmp <- rep(c(4:c(8)), each=4)
  tmp[seq(1, 17, 4)] <- 0
  lmat <- matrix(c(0,2,2,2,3,1,1,1, tmp, 0,9,10,11), byrow=F, nrow=4)
  lwid <- c(0.5, 3, rep(0.05, 5), 0.5)  # values of widths of columns in the layout
  lhei <- c(0.5, 1, 1, 1)     # values of heigts of rows on the device
  # set such that first row and column (dendro) is half size of rest
  # there's 2 narrow columns more than there are monophyly statistics, because
  # those are sandwiched between two white columns
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  
  # plot parameters
  heatMargins <- c(0, 1)
  phyloMargins <- c(2, 1, 2, 1)
  labCol <- NULL
  labRow <- NULL
  cexCol <- 0.35
  cexRow <- 0.35
  
  # color parameters
  # for heatmap
  tmp <- c("white", brewer.pal(9, "YlOrRd"))    ## for now, from yellow to red via orange
  colsPalette <- colorRampPalette(tmp, interpolate="linear", bias=0.4)
  cols <- (colsPalette(100))
  
  # data properties
  x <- mat
  di <- dim(x)
  nr <- di[1L]
  nc <- di[2L]

  ## get the dendrograms and reordering indices
   hcr <- hclust(dist(x, method="euclidean"), method="average") # average for a UPGMA tree
   ddr <- as.dendrogram(hcr)
   ddr <- dendsort::dendsort(ddr)
   ddc <- ddr
   rowInd <- order.dendrogram(ddr)
   colInd <- rowInd
  # 
  # ## reorder data 
   x <- x[rowInd, colInd] 
 
  #  reorder the rows of the image matrix by pp for AC and plot
  # pps_AC <- pps[rownames(mat), "pp_AC"]
  # pps_AP <- pps[rownames(mat), "pp_AP"]
  # pps_PC <- pps[rownames(mat), "pp_PC"]
  # 
  # p <- data.frame(pps_AC=pps_AC, pps_AP=pps_AP, pps_PC=pps_PC)
  # 
  # # re-order
  # strategy: 
  # 1. find the best-supported topology for each gene
  # 2. for each gene in turn, find the order of support values
  # 3. stack the rows by genesupport then 
  
  # helper function to get the max number for all rows of a dataframe, giving NA when its all-tied.
  # rowWhich.max <- function(dataframe) {
  #   result <- NULL
  #   for(i in 1:nrow(dataframe)){
  #     row <- unname(dataframe[i, ])
  #     if(length(unique(c(row)))  != 1 ) { # there's no competing "best" values
  #       thisRep <- which.max(row)
  #     } else { thisRep <- NA} # there's equal support all around
  #     result <- c(result, thisRep)
  #   }
  #   return(result)
  # }
  # # use helperfunction to get the optimal reordering sequence
  # whichtree <- rowWhich.max(p)
  # p <- cbind(p, whichtree)
  # p1 <- p[which(p$whichtree == 1), c(1:3)]
  # p2 <- p[which(p$whichtree == 2), c(1:3)]
  # p3 <- p[which(p$whichtree == 3), c(1:3)]
  # pNA<- p[which(is.na(p$whichtree)), c(1:3)]
  # 
  # p1 <- p1[order(p1[,1]), ]
  # p2 <- p2[order(p2[,2]), ]
  # p3 <- p3[order(p3[,3]), ]
  # 
  # reorderedPP <- rbind(p1, p2, p3, pNA)
  # 
  # rI <- reorderIndex <- as.numeric(rownames(reorderedPP))
  # # length(rI) == length(unique(rI)) # worked
  # # length(rI) == nrow(mat) # worked
  # 
  # # use reorderIndex to reorder mat
  # x <- mat[rI, rI]
  # di <- dim(x)
  # nr <- di[1L]
  # nc <- di[2L]
  
  ## 1. plot the main bit
  par(mar = c(heatMargins[1L], 0, 0, heatMargins[2L]))
  image(1L:nc, 1L:nr, x, xlim = 0.5+ c(0, nc), ylim = 0.5+ c(0, nr),
        axes = FALSE, xlab = "", ylab = "", col=cols)
  
  # 1b. add locus names
  
  rowlabs <- rownames(mat)[rowInd] # names of loci ordered from bottom to top 
  # get their position in allLocusResults
  # for each element of rowlabs, find the corresponding position in allLocusResults
  position <- NULL
  for(i in rowlabs) {position <- c(position, which(rownames(allLocusResults) == i))}
  rowlabs <- sub("c2p1pgts2_", "", rowlabs)
  rowlabs <- sub("_mb", "", rowlabs)
  # combine locus number (as in the tree plots) and the locus name
  rowlabs <- paste(rowlabs, "(#", sprintf("%3i", position), ")", sep="")
  # use this as an axis label  
  
  axis(2, 1L:nc, labels = rowlabs, las = 2, line = -1, tick = 0,   ## to plot labels
       cex.axis = 0.4, family="mono")
  
  ## 2, nothing
  frame()  # skip this plot; space was used for the axis labels
  
  #  3. the dendrogram :
  par(mar = c(0, 0, 0, heatMargins[2L]))
  plot(ddr,		axes = FALSE, xaxs = "i", leaflab = "none")
  
  ## the support symbols
  par(mar = c(0, 0.1, 0, 0.1))
  
  # the relevant locus results, indexed by heatmaprows
  n <- nrow(x)
  locNameInd <- rownames(mat)[rowInd] ## order of loci from bottom to top 
  
  pp_cols <- brewer.pal(9, "Greys")[2:8]
  pp_cols<- c(rep(pp_cols, each=2), "black")  # 15 indices of near white to black; 
  
  # first a totally white column
  image(rbind(1L:n), col = "white", axes = FALSE)
  # then a column for support values;
  # get pps from allLocusResults, index via locusnames.
  AP_pp <- allLocusResults[locNameInd, "pp_AP"]
  image(rbind(1L:n), col = pp_cols[(AP_pp*14)+1], axes = FALSE)
  mtext("AP", las=2, cex=0.6)
  AC_pp <- allLocusResults[locNameInd, "pp_AC"]   
  image(rbind(1L:n), col = pp_cols[(AC_pp*14)+1], axes = FALSE)
  mtext("AC", las=2, cex=0.6)
  PC_pp <- allLocusResults[locNameInd, "pp_PC"]
  image(rbind(1L:n), col = pp_cols[(PC_pp*14)+1], axes = FALSE)
  mtext("PC", las=2, cex=0.6)
  # then another totally white column
  image(rbind(1L:n), col = "white", axes = FALSE)
  
  ## finally, plot the trees
  
  par(mar = c(2,0,2,0))
  # get the trees to draw:
  APtree <- read.tree(text="(Cact.,(Anac.,Port.)0.519);")
  ACtree <- read.tree(text="(Port.,(Anac.,Cact.)0.251);")
  PCtree <- read.tree(text="(Anac.,(Port.,Cact.)0.230);")
  
  plot(ACtree, label.offset=0.2, edge.width=3, cex=1.2)
  drawSupportOnEdges(ACtree$node.label, adj=c(0.5, -0.5), frame="none")
  
  plot(PCtree, label.offset=0.2, edge.width=3, cex=1.2)
  drawSupportOnEdges(PCtree$node.label, adj=c(0.5, -0.5), frame="none")
  
  plot(APtree, label.offset=0.2, edge.width=3, cex=1.2)
  drawSupportOnEdges(APtree$node.label, adj=c(0.5, -0.5), frame="none")
  
  dev.off()
  
  
  # all done.  Yay!
 
  