mytree.asymmetric.age <-
function (age, waitsp, waitext, complete=TRUE, tiplabel, shiftsp, shiftext, sampling)
{ 
  # create waiting time functions #
  if (is.function(waitsp)){
    rnumbsp <- waitsp#parse(text="waitsp()")
  } else if (is.character(waitsp)){
    rnumbsp <- express.distribution(waitsp)
  }
  
  if (is.function(waitext)){
    rnumbext <- waitext#parse(text="waitext()")
    firstextpar <- "funk"
  } else if (is.character(waitext)){
    rnumbext <- express.distribution(waitext)
    firstextpar <- get.first.par.distribution(waitext)
  }
  
  if (is.function(shiftsp$strength)){
    rnumbshiftsp <- shiftsp$strength#parse(text="shiftsp$strength()")
  } else if (is.character(shiftsp$strength)){
    rnumbshiftsp <- express.distribution(shiftsp$strength)
  }
  
  if (is.function(shiftext$strength)){
    rnumbshiftext <- shiftext$strength#parse(text="shiftext$strength()")
  } else if (is.character(shiftext$strength)){
    rnumbshiftext <- express.distribution(shiftext$strength)
  }
  
  # end create waiting time functions #
  
  
  labellivingsp=tiplabel[1]
  labelextinctsp=tiplabel[2]
  shiftsplabel=tiplabel[3]
  shiftextlabel=tiplabel[4]
  
  shiftspprob=shiftsp$prob
  shiftextprob=shiftext$prob 
  
  #fucntion for testing for shift speciation
  testshiftsp <- function (spt){
    if (shiftspprob!=0){
      if (runif(1,0,1)<shiftspprob) { #cheking if shift happens based on a uniform distribution U[0,1]
        shiftspstrength <- rnumbshiftsp()#eval(rnumbshiftsp) #draw shift strength
        #print(paste("Hey, a speciation shift happened at node", nextsp - i, "with strenght", shiftspstrength ))
      }else{
        shiftspstrength <- shiftspm[shiftspm[,"node"]==species,"strength"] #inherit ancestor shift
        #print(paste("No speciation shift!",nextsp - i, "is inheriting the shift of",shiftspstrength, "from", species ))
      }    
    }else{
      shiftspstrength <- 1
    }
    spt <- spt*shiftspstrength #update speciation time
    shiftspmr <- c((nextsp - i),shiftspstrength ) #link shift with node to be added to shiftspm...
    attributes(spt) <- NULL
    attributes(shiftspmr) <- NULL
    return(list(spt=spt,shiftspmr=shiftspmr))#return updated shift and updated shiftspmatrix
  }
  
  #fucntion for testing for shifts on extinction process
  testshiftext <- function (extt){
    if (shiftextprob!=0){
      if (runif(1,0,1)<shiftextprob) { #cheking if shift happens based on a uniform distribution U[0,1]
        shiftextstrength <- rnumbshiftext()#eval(rnumbshiftext) #draw shift strength
        #print(paste("Hey, an extinction shift happened at node", nextsp - i, "with strenght", shiftextstrength ))
      }else{
        shiftextstrength <- shiftextm[shiftextm[,"node"]==species,"strength"] #inherit ancestor shift
        #print(paste("No extinction shift!",nextsp - i, "is inheriting the shift of",shiftextstrength, "from", species ))
      }    
    }else{
      shiftextstrength <- 1
    }
    extt <- extt*shiftextstrength #update extinction time
    shiftextmr <- c((nextsp - i),shiftextstrength ) #link shift with node to be added to shiftextm...
    attributes(extt) <- NULL
    attributes(shiftextmr) <- NULL
    return(list(extt=extt,shiftextmr=shiftextmr))#return updated shift and updated shiftextmatrix
  }
  ##SM-E
  
  trajectory  <- function (trace){
    #trace should indicate the edge number that is to be followed until the origin
    trajectory  <- NULL
    while ( length( which(edge[,2] == trace)) ){
      atual  <- which(edge[,2] == trace)
      trajectory  <- c(edge.length[atual], trajectory)
      trace  <- edge[atual,1]
    }
    return(trajectory)
  }  
  # tracing back the time for one species until the "until"
  trajectoryuntil  <- function (trace, until){
    #trace should indicate the edge number that is to be followed until the "until"
    atual <- trace
    trajectoryuntil  <- NULL
    while ( trace != until ){
      atual  <- which(edge[,2] == trace)
      trajectoryuntil  <- c(edge.length[atual], trajectoryuntil)
      trace  <- edge[atual,1]
    }
    return(trajectoryuntil)
  }
  

  stop <- FALSE
  mytree <-list(edge=NULL, tip.label=NULL, edge.length=NULL, Nnode=NULL, root.edge=NULL, age=NULL, 
                shiftsp=NULL,shiftext=NULL, shifted.sp.living=NULL,shifted.sp.extinct=NULL,
                shifted.ext.living=NULL, shifted.ext.extinct=NULL)
  class(mytree) <- "phylo"
  edge <- matrix(c(-1,-2), ncol=2)
  leaves <- NULL
  realleaves <- NULL
  extinct <- NULL
  tip.label <- NULL
  ## SM - S
  shiftedspliving <- NULL
  shiftedextliving <- NULL
  shiftedspextinct <- NULL
  shiftedextextinct <- NULL
  ## SM - E

  
  ##SM-S creating shiftsp and shiftext matrixes...
  shiftspm <- matrix(c(-1,1,-2,1), byrow=TRUE, ncol=2, dimnames=list(NULL,c("node","strength")))
  shiftextm <- matrix(c(-1,1,-2,1), byrow=TRUE, ncol=2, dimnames=list(NULL,c("node","strength")))
  ##SM-E
  
  # initial if in case the (-1,-2) edge get extinct or bigger than age
  spt <- rnumbsp()#eval(rnumbsp)
  #to remove NaN warnings messages in case of zero inside the distribution parameters
  {
  if (firstextpar == 0)
  {
    extt <- suppressWarnings(rnumbext())#eval(rnumbext)
  }
  else
  {
    extt <- rnumbext()#eval(rnumbext)
  }
  }
  {
  if (is.nan(extt))
  {
    extt <- age + 0.1 # we add 0.1 to age so that an extinction event will never happen
  }
  }
  {
  if (spt <= extt)
  {
    status <- "sp" #occurred an speciation
    edge.length <- spt
    leaves <- -2
  }
  else
  {
    status <- "ext" #occurred an extinction
    edge.length <- extt
    extinct <- -2
    stop <- TRUE
  }
  }
  
  {
  if (min(spt,extt) >= age)
  {
    edge.length <- age
    stop <- TRUE
    {
    if (status == "sp")
    {
      realleaves <- leaves
      leaves <- NULL					
    }
    else
    {
      realleaves <- extinct			
      extinct <- NULL
    }
    }
  }
  }
  # in here are stored all ages for death of species and their descendent and status
  godbook <- list(from=-1,to=-2,destiny=extt, status=status)
  
  # for the while (increase of the tree)
  while (stop == FALSE) 
  {
    species <- leaves[1]
    nextsp <- min(edge[,2])
    #print(paste("working on descendantes from", species))
    
    #for the mother prolong
    i <- 1
    edge <- rbind( edge, c(species, (nextsp - i)))
    motherage <- sum(    trajectoryuntil(   godbook$to[which(godbook$to == species)]  ,    godbook$from[which(godbook$to == species)]))
    spt <- rnumbsp()#eval(rnumbsp)
    #SM-S
    testshiftspout <- testshiftsp(spt)
    spt <- testshiftspout$spt
    ##SM -E
    # here we sample until the speciation (i.e. spt) is bigger than the traject length of the leaf, because the process should be valid for speciation as for extinction.
    #
    ##
    #probablydelta=FALSE
    count <- 0
    while (spt <= motherage) {
      # print(paste("again.. because spt was=", spt))
      spt <- rnumbsp()#eval(rnumbsp)
      #SM-S
      testshiftspout <- testshiftsp(spt)
      spt <- testshiftspout$spt
      count <- count+1
      if (count>10 & spt==motherage) {
        spt <- rnumbext() +1
        break()
      }
      ##SM -E
    }
    ## SM - S
    shiftspm <- rbind(shiftspm, testshiftspout$shiftspmr)
    #updating extinction strenght.. inheriting from previous node = "mother"
    shiftextm <- rbind(shiftextm ,c((nextsp - i),(shiftextm[shiftextm[,"node"]==species,"strength"]))) #inherit ancestor shift)
    ##SM - E
    # here we cut the overal sampled spt to add only the additional age
    spt <- spt-motherage
    ##
    #
    #on the case of the mother, we don't sample any more extinction time, since it has been given already
    {
    if (motherage + spt <= godbook$destiny[which(godbook$to == species)])
    {
      status <- "sp" #occurred an speciation
      edge.length <- c(edge.length, spt)
      leaves <- c(leaves, (nextsp - i))
      godbook$to[which(godbook$to == species)] <- nextsp -i
    }
    else
    {
      status <- "ext" #occurred an extinction
      edge.length <- c(edge.length, godbook$destiny[which(godbook$to == species)]  - sum(  trajectoryuntil(species,  godbook$from[which(godbook$to == species)])))
      extinct <- c(extinct, (nextsp - i))
      godbook$status[which(godbook$to == species)] <- status
    }
    }
    #
    traject <- trajectory(nextsp - i)
    {
    if ( sum(traject) >= age) 
    {
      edge.length[length(edge.length)] <- age - sum(traject[1: (length(traject)-1) ] ) # reduce the last edge.length to make the sum be equals age
      # now we see if it extincted or speciated to know from were to take it out and put into the realleaves
      {
      if (status == "sp")
      {
        realleaves <- c( realleaves, leaves[length(leaves)])
        leaves <- leaves[-length(leaves)]					
      }
      else
      {
        realleaves <- c( realleaves, extinct[length(extinct)])			
        extinct <- extinct[-length(extinct)]
      }
      }
    }  
    }
    #for the child birth
    i <- 2
    edge <- rbind( edge, c(species, (nextsp - i)))
    spt <- rnumbsp()#eval(rnumbsp)
    ##SM -S now we examine for shifts and update spt and shiftspm
    testshiftspout <- testshiftsp(spt)
    spt <- testshiftspout$spt
    shiftspm <- rbind(shiftspm, testshiftspout$shiftspmr)
    ##SM -E
    #to remove NaN warnings messages
    {
    if (firstextpar == 0)
    {
      extt <- suppressWarnings(rnumbext())#eval(rnumbext)
    }
    else
    {
      extt <- rnumbext()#eval(rnumbext)
    }
    }
    # if to see if the user simulates with extinction = ZERO and avoid error generated by expression-distribution when rate equals zero
    {
    if (is.nan(extt))
    {
      extt <- age +1 #  we add 1 to age so that an extinction event will never happen
    }
    }
    
    ##SM -S now we examine for shifts and update extt and shiftextm
    testshiftextout <- testshiftext(extt)
    extt <- testshiftextout$extt
    shiftextm <- rbind(shiftextm, testshiftextout$shiftextmr)
    ##SM -E
    
    godbook$from <- c(godbook$from, species)
    godbook$to <- c(godbook$to, nextsp - i)	
    godbook$destiny <- c(godbook$destiny, extt)
    {
    if (spt <= extt)
    {
      status <- "sp" #occurred an speciation
      edge.length <- c(edge.length, spt)
      leaves <- c(leaves, (nextsp - i))
    }
    else
    {
      status <- "ext" #occurred an extinction
      edge.length <- c(edge.length, extt)
      extinct <- c(extinct, (nextsp - i))
    }
    }
    godbook$status <- c(godbook$status, status)
    #
    traject <- trajectory(nextsp - i)
    {
    if ( sum(traject) >= age) 
    {
      edge.length[length(edge.length)] <- age - sum(traject[1: (length(traject)-1)]) # reduce the last edge.length to make the sum be equals age
      # now we see if it extincted or speciated to know from were to take it out and put into the realleaves
      {
      if (status == "sp")
      {
        realleaves <- c( realleaves, leaves[length(leaves)])
        leaves <- leaves[-length(leaves)]					
      }
      else
      {
        realleaves <- c( realleaves, extinct[length(extinct)])			
        extinct <- extinct[-length(extinct)]
      }
      }
    }  
    }	
    
    leaves <- leaves[-1]
    {
    if (length(leaves) == 0)
    {
      stop <- TRUE	
    }
    }
  }
  
  
  
  
  # final if... in case of (stop == TRUE) , we write the tree "mytree"
  {
  if (stop == TRUE) 
  {
    #### replacing to the ape format
    prealleaves <- realleaves
    {
    if (length(realleaves) > 0)
    {
      realleaves <- c(1:length(realleaves))
      i <- 1
      for (i in 1:length(realleaves))
      {
        edge[ which(edge[,2] == prealleaves[i]), 2 ] <- realleaves[i]
        ## SM- S changing nodes names
        shiftspm[ which(shiftspm[,"node"] == prealleaves[i]), "node" ] <- realleaves[i]
        shiftextm[ which(shiftextm[,"node"] == prealleaves[i]), "node" ] <- realleaves[i]
        ## SM- E
      }
      tip.label <- paste(labellivingsp, realleaves, sep = "")
      ## SM - S creating vector with shifted species and updating shifted names
      #for shifts on speciation
      shiftedspliving <- rep(0,length(realleaves)) #preparing vector
      
      #shiftedspliving[shiftspm[shiftspm[,"node"]%in%realleaves,"strength"]!=1] <- 1
      
      data <- shiftspm[shiftspm[,"node"]%in%realleaves,]

      if (length(shiftedspliving)==1){
        data <- t(as.matrix(data))
      }
      shiftedspliving[data[sort.list(data[,"node"]),"strength" ]!=1] <- 1
      
      tip.label[shiftedspliving==1] <- paste(tip.label[shiftedspliving==1], shiftsplabel, sep=" ")
      shiftedspliving <- cbind(realleaves,shiftedspliving)
      colnames(shiftedspliving) <- c("LivingSpecies", "shift")
      #for shifts on extinction
      shiftedextliving <- rep(0,length(realleaves))
      
      #shiftedextliving[shiftextm[shiftextm[,"node"]%in%realleaves,"strength"]!=1] <- 1
      
      data <- shiftextm[shiftextm[,"node"]%in%realleaves,]

      if (!is.matrix(data)){
        data <- t(as.matrix(data))
      }
      shiftedextliving[data[sort.list(data[,"node"]),"strength" ]!=1] <- 1
      
      tip.label[shiftedextliving==1] <- paste(tip.label[shiftedextliving==1], shiftextlabel, sep=" ")
      shiftedextliving <- cbind(realleaves,shiftedextliving)
      colnames(shiftedextliving) <- c("LivingSpecies", "shift")
      ## SM - E
    }
    }
    pextinct <- extinct
    {
    if (length(extinct) > 0)
    {
      extinct <- c((length(realleaves)+1):(length(realleaves)+length(extinct)))
      i <- 1
      for (i in 1:length(extinct))
      {
        edge[ which(edge[,2] == pextinct[i]), 2 ] <- extinct[i]
        ## SM- S changing nodes names
        shiftspm[ which(shiftspm[,"node"] == pextinct[i]), "node" ] <- extinct[i]
        shiftextm[ which(shiftextm[,"node"] == pextinct[i]), "node" ] <- extinct[i]
        ## SM- E
      }
      tip.label.tail <- paste(labelextinctsp, extinct, sep = "")
      ## SM - S creating vector with shifted species and updating shifted names
      #for shifts on speciation
      shiftedspextinct <- rep(0,length(extinct))
      data <- shiftspm[shiftspm[,"node"]%in%extinct,]
      #error handling in case only one species got extinct and data is transformed into a vector and not a matrix...
      if(class(data)=="numeric"){
        data <- matrix(data, ncol=2)
      }
      shiftedspextinct[data[sort.list(data[,1]),2 ]!=1] <- 1
      tip.label.tail[shiftedspextinct==1] <- paste(tip.label.tail[shiftedspextinct==1], shiftsplabel, sep=" ")
      shiftedspextinct <- cbind(extinct,shiftedspextinct)
      colnames(shiftedspextinct) <- c("ExtinctSpecies", "shift")
      #for shifts on extinction
      shiftedextextinct <- rep(0,length(extinct))
      data <- shiftextm[shiftextm[,"node"]%in%extinct,]
      #error handling in case only one species got extinct and data is transformed into a vector and not a matrix...
      if(class(data)=="numeric"){
        data <- matrix(data, ncol=2)
      }
      shiftedextextinct[data[sort.list(data[,1]),2 ]!=1] <- 1
      tip.label.tail[shiftedextextinct==1] <- paste(tip.label.tail[shiftedextextinct==1], shiftextlabel, sep=" ")
      shiftedextextinct <- cbind(extinct,shiftedextextinct)
      colnames(shiftedextextinct) <- c("ExtinctSpecies", "shift")
      ## SM - E
      tip.label <- c(tip.label, tip.label.tail )
      
    }
    }
    #regarding the edges that lead to an extinct or leaving final species, but are not the final edges
    potheredges <- levels(as.factor(edge[edge <0]))
    otheredges <- rev(seq((max(realleaves, extinct)+1), length.out=length(potheredges)))
    i <- 1
    for (i in 1:length(potheredges))
    {
      edge[ edge == potheredges[i] ] <- otheredges[i]
      ## SM- S changing nodes names
      shiftspm[shiftspm[,"node"]==potheredges[i],"node"] <- otheredges[i]
      shiftextm[shiftextm[,"node"]==potheredges[i],"node"] <- otheredges[i]
      ## SM- E
    }
    mytree$edge <- edge
    mytree$tip.label <- tip.label
    mytree$edge.length <- edge.length
    mytree$Nnode <-  length(realleaves) + length(extinct)
    mytree$root.edge <- edge.length[1]
    mytree$age <- age
    ## SM - S
    mytree$shiftsp <- shiftspm
    mytree$shiftext <- shiftextm
    mytree$shifted.sp.living <- shiftedspliving
    mytree$shifted.sp.extinct <- shiftedspextinct
    mytree$shifted.ext.living <- shiftedextliving
    mytree$shifted.ext.extinct <- shiftedextextinct
    ## SM - E
  }	
  }
  
  #final handling before plotting to handle ape limitations
  {
  if (length(realleaves) == 0)
  {
    # in case no specie is surviving until final simulation time
    mytree <- 0
  }
  else
  {
    {
      if ( length(realleaves) == 1 & complete == FALSE)
      {
        #in case only one specie is surviving, even if other speciations events occurred in the history
        mytree <- 1
      }
      else
      {
        {
          if (length(realleaves)==1 & length(extinct)==0 & complete==TRUE)
          {
            #in case only one species is surviving and was the only one that existed
            mytree <- 1
          }
          else
          {
            #in case non of the above condition is fulfilled, there will be a tree with no initial branch, tree starts at the MRCA
            #this is done to be aple to plot with the`ape` package 
            mytree <- collapse.singles(mytree)
            #checking status of 'complete' and take or dont take extincted species out of final tree
            {
            if (complete == FALSE)	
            {
              mytree<- drop.fossil(mytree)
              mytree$root.edge <- age - max(getx(mytree)) # restore back the root length!
            }
            }		
          }
        }		
      }
    }	
  }
  }
  if (sampling$frac!=1 & class(mytree)=="phylo"){ # do sampling.....
    mytree <- sample.mytree(mytree.=mytree, realleaves.=realleaves, extinct.=extinct, sampling.=sampling)
  }
  return(mytree)
}
