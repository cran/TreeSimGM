track.shift <- function(shift="sp",tree, node){ #track shift="sp" or "ext"
  if (!is.numeric(node)){
    nodeold <- node
    node <- as.numeric(gsub("[^0-9]","",nodeold))
    print(paste0("we identified that '",nodeold, "' refers to node number ", node ))
  }
  edge <- tree$edge
  shiftspm <- tree$shiftsp
  shiftextm <- tree$shiftext
  
  ###################################################### inner functions ######################################
  # tracing back the shifts on speciation for one species until origin, name of the attributes tell the node which change happened
  trajectoryshiftsp  <- function (trace){
    #trace should indicate the edge number that is to be followed until the origin
    trajectory  <- NULL
    while ( length( which(edge[,2] == trace)) ){
      atual  <- which(edge[,2] == trace)
      add <- shiftspm[shiftspm[,"node"]==edge[atual,2],"strength"]
      attributes(add)$names <- shiftspm[shiftspm[,"node"]==edge[atual,2],"node"]
      trajectory  <- c(add, trajectory)
      trace  <- edge[atual,1]
    }
    return(trajectory)
  }
  # tracing back the shifts on extinction for one species until origin, name of the attributes tell the node which change happened
  trajectoryshiftext  <- function (trace){
    #trace should indicate the edge number that is to be followed until the origin
    trajectory  <- NULL
    while ( length( which(edge[,2] == trace)) ){
      atual  <- which(edge[,2] == trace)
      add <- shiftextm[shiftextm[,"node"]==edge[atual,2],"strength"]
      attributes(add)$names <- shiftextm[shiftextm[,"node"]==edge[atual,2],"node"]
      trajectory  <- c(add, trajectory)
      trace  <- edge[atual,1]
    }
    return(trajectory)
  }
  ############################################ end inner functions ######################################
  
  if (shift=="sp"){
    print("sp")
    trajectory <- trajectoryshiftsp(trace=node)
  }else if (shift=="ext"){
    print("ext")
    trajectory <- trajectoryshiftext(trace=node)
  }else{
    print("please choose either `sp` (tracking shifts on speciation) or `ext` (tracking shifts on extinction) for the shift parameter")
  }
  
  return(trajectory)
  
}


######### Support Functions ##############
express.distribution <- function(distributionstring){
  if (!is.character(distributionstring)){
    stop("please provide a string as waiting time or shift strength")
  }
  split <- strsplit(distributionstring, "(" , fixed  =T)
  return(parse(text=paste0(split[[1]], collapse = "(1,")))
}

get.first.par.distribution <- function(distributionstring){
  split <- strsplit(distributionstring, "(" , fixed  =T)
  split1 <- strsplit(split[[1]][2], "," , fixed  =T)
  splot2 <- as.numeric(gsub(')|,','',split1[[1]][1]))
  return(splot2)
}

ape.mod.drop.tip <- function (phy, tip,  root.edge = 0)
  #adapted from ape citation("ape")
  # Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of phylogenetics and evolution
  # in R language. Bioinformatics 20: 289-290.
{
  Ntip <- length(phy$tip.label)
  if (is.character(tip)) 
    tip <- which(phy$tip.label %in% tip)
  
  out.of.range <- tip > Ntip
  if (any(out.of.range)) {
    # warning("some tip numbers were larger than the number of tips: they were ignored")
    tip <- tip[!out.of.range]
  }
  if (!length(tip)) 
    return(phy)
  if (length(tip) == Ntip) {
    #warning("drop all tips of the tree: returning NULL")
    return(0)
  }
  wbl <- !is.null(phy$edge.length)
  if (length(tip) == Ntip - 1 ) {
    i <- which(phy$edge[, 2] == (1:Ntip)[-tip])
    res <- list(edge = matrix(2:1, 1, 2), tip.label = phy$tip.label[phy$edge[i, 
                                                                             2]], Nnode = 1L)
    class(res) <- "phylo"
    if (wbl) 
      res$edge.length <- phy$edge.length[i]
    if (!is.null(phy$node.label)) 
      res$node.label <- phy$node.label[phy$edge[i, 1] - 
                                         Ntip]
    return(res)
  }
  
  phy <- reorder.phylo(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- dim(phy$edge)[1]
  
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !logical(Nedge)
  keep[match(tip, edge2)] <- FALSE
  
  ints <- edge2 > Ntip
  repeat {
    sel <- !(edge2 %in% edge1[keep]) & ints & keep
    if (!sum(sel)) 
      break
    keep[sel] <- FALSE
  }
  
  if (root.edge && wbl) {
    degree <- tabulate(edge1[keep])
    if (degree[ROOT] == 1) {
      j <- integer(0)
      repeat {
        i <- which(edge1 == NEWROOT & keep)
        j <- c(i, j)
        NEWROOT <- edge2[i]
        degree <- tabulate(edge1[keep])
        if (degree[NEWROOT] > 1) 
          break
      }
      keep[j] <- FALSE
      if (length(j) > root.edge) 
        j <- 1:root.edge
      NewRootEdge <- sum(phy$edge.length[j])
      if (length(j) < root.edge && !is.null(phy$root.edge)) 
        NewRootEdge <- NewRootEdge + phy$root.edge
      phy$root.edge <- NewRootEdge
    }
  }
  
  if (!root.edge) 
    phy$root.edge <- NULL
  phy$edge <- phy$edge[keep, ]
  if (wbl) 
    phy$edge.length <- phy$edge.length[keep]
  TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
  oldNo.ofNewTips <- phy$edge[TERMS, 2]
  
  n <- length(oldNo.ofNewTips)
  phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
  phy$tip.label <- phy$tip.label[-tip]
  
  phy$Nnode <- dim(phy$edge)[1] - n + 1L
  newNb <- integer(Ntip + Nnode)
  newNb[NEWROOT] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
  phy$edge[, 1] <- newNb[phy$edge[, 1]]
  storage.mode(phy$edge) <- "integer"
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
  collapse.singles(phy)
}


sample.mytree <- function(mytree., realleaves., extinct., sampling., complete.) {
  if (sampling.$frac<0 | sampling.$frac>1) {
    stop("sampling fraction needs to be between 0 and 1")
  }
  # frac=0.8 # sampling fraction!
  alltips <- length(mytree.$tip.label)
  
  # howmanytips will be removed
  ndelete <- round((1 - sampling.$frac) * alltips)
  
  if (sampling.$branchprop==TRUE) {
    alltips.branches <- mytree.$edge[,2]<=alltips
    alltips.numbers <- mytree.$edge[alltips.branches,2]
    alltips.branches <- mytree.$edge.length[alltips.branches]
    names(alltips.branches) <- alltips.numbers
    alltips.branches <- alltips.branches/sum(alltips.branches)
    alltips.branches <- 1-alltips.branches
    alltips.branches <- alltips.branches/sum(alltips.branches)
    delete <- as.integer(names(sample(x=alltips.branches, size=ndelete, prob=alltips.branches)))
  } else {
    delete <- sample(x=mytree.$tip.label, size=ndelete)
  }
  # do tip_dropping...adapted from ape
  sampledtree <- ape.mod.drop.tip(mytree., tip=delete)
  
  if (!is.numeric(sampledtree)) {
    tipsfromlabel <- regmatches(sampledtree$tip.label, gregexpr("[[:digit:]]+", sampledtree$tip.label))
    tipsfromlabel <- as.numeric(unlist(tipsfromlabel))
    Nliving <- realleaves.%in%tipsfromlabel
    if ( sum(Nliving)==0 ) {
      sampledtree <- 0
    } else if ( sum(Nliving)==1 & complete.==TRUE & sum(extinct.%in%tipsfromlabel)==0 ) {
      sampledtree <- 1
    }
  }
  # clean sampled tree...
  if (!is.numeric(sampledtree)) {
    sampledtree$shiftsp <- NULL
    sampledtree$shiftext <- NULL
    sampledtree$shifted.sp.living <- NULL
    sampledtree$shifted.sp.extinct <- NULL
    sampledtree$shifted.ext.living <- NULL
    sampledtree$shifted.ext.extinct <- NULL
  }
  suppressWarnings(sampledtree$beforesampling <- mytree.)
  mytree. <- sampledtree
  return(mytree.)
}

