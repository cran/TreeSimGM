mytree.asymmetric.taxa <-
function (m, distributionspname, distributionspparameters, distributionextname, distributionextparameters, complete=TRUE, labellivingsp="sp.", labelextinctsp="ext.", 
          shiftspprob, shiftdistributionspname, shiftdistributionspparameters, shiftextprob, shiftdistributionextname, shiftdistributionextparameters, shiftsplabel, shiftextlabel)
{ 

udistributionspparameters <- capture.output (cat(distributionspparameters, sep=","))
udistributionextparameters <- capture.output (cat(distributionextparameters, sep=","))
rnumbsp <- parse(text=paste(distributionspname, "(1,", udistributionspparameters,")"))
rnumbext <- parse(text=paste(distributionextname, "(1,", udistributionextparameters,")"))
# for a delta test or dirac delta function test of the algorithim
#rnumbsp <- 0.5
#rnumbext <- 0.6

##SM-S
ushiftdistributionspparameters <- capture.output (cat(shiftdistributionspparameters, sep=","))
ushiftdistributionextparameters <- capture.output (cat(shiftdistributionextparameters, sep=","))
rnumbshiftsp <- parse(text=paste(shiftdistributionspname, "(1,", ushiftdistributionspparameters,")"))
rnumbshiftext <- parse(text=paste(shiftdistributionextname, "(1,", ushiftdistributionextparameters,")"))

#testing for shift speciation
testshiftsp <- function (spt){
  if (shiftspprob!=0){
    if (runif(1,0,1)<shiftspprob) { #cheking if shift happens based on a uniform distribution U[0,1]
      shiftspstrength <- eval(rnumbshiftsp) #draw shift strength
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


#testing for shifts on extinction process
testshiftext <- function (extt){
  if (shiftextprob!=0){
    if (runif(1,0,1)<shiftextprob) { #cheking if shift happens based on a uniform distribution U[0,1]
      shiftextstrength <- eval(rnumbshiftext) #draw shift strength
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

# tracing back the time for one species until origin
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
bingo <- FALSE

while (bingo == FALSE)	# while bingo start
{
  ##SM - S creating shiftsp and shiftext matrixes...
  shiftspm <- matrix(c(-1,1,-2,1), byrow=TRUE, ncol=2, dimnames=list(NULL,c("node","strength")))
  shiftextm <- matrix(c(-1,1,-2,1), byrow=TRUE, ncol=2, dimnames=list(NULL,c("node","strength")))
  shiftedspliving <- NULL
  shiftedspextinct <- NULL
  shiftedextliving <- NULL
  shiftedextextinct <- NULL
  ##SM - E
  stop <- FALSE
	stopsearch <- FALSE
	mytree <-list(edge=NULL, tip.label=NULL, edge.length=NULL, Nnode=NULL, root.edge=NULL, age=NULL, 
	              shiftsp=NULL,shiftext=NULL, shifted.sp.living=NULL,shifted.sp.extinct=NULL,
	              shifted.ext.living=NULL, shifted.ext.extinct=NULL)
	class(mytree) <- "phylo"
	edge <- matrix(c(-1,-2), ncol=2)
	leaves <- NULL
	realleaves <- NULL
	extinct <- NULL
	realextinct <- NULL
	age <- NULL
	timeline <- NULL
	extinct <- NULL
	tip.label <- NULL
	pnodges <- NULL
	timeline <- NULL

# initial if in case the (-1,-2) edge get extinct or bigger than age
spt <- eval(rnumbsp)
#to remove NaN warnings messages
{
if (distributionextparameters[1] == 0)
{
	extt <- suppressWarnings(eval(rnumbext))
}
else
{
	extt <- eval(rnumbext)
}
}
# if to see if the user simulates with extinction = ZERO and avoid error generated by exp. distribution when rate equals zero
{
if (is.nan(extt))
{
	extt <- NA # we add so that the sp will always occurs, and we will never have an extinction (later on, NA will be replaced by the age)
}
}

{
if (spt <= extt | is.na(extt))
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

godbook <- list(from=-1,to=-2, destiny=extt, age=edge.length[1], status=status)
pnodges <- c(leaves, extinct)
timeline <-edge.length[1]

# for the while (increase of the tree)
while (stop == FALSE) 
{
	stopsearch <- FALSE
	positionnextsp <- which(timeline == min(timeline[pnodges%in%leaves]))[1] #which(timeline == min(timeline))[1] #ATTENTION, we take the first element since it can happen that we get two equal numbers
	nextsp <- min(edge[,2])
	species <- pnodges[positionnextsp] #print(paste("Working on descendents from species",species))
	ptime <- timeline[positionnextsp] #ptime stands for PreviousTIME
	pnodges <- pnodges[-positionnextsp]
	timeline <- timeline[-positionnextsp]
	
	#for the mother to prolong
	i <- 1
	edge <- rbind( edge, c(species, (nextsp - i)))
	spt <- eval(rnumbsp)
	##SM -S now we examine for shifts and update spt and shiftspm
	testshiftspout <- testshiftsp(spt)
	spt <- testshiftspout$spt
	##SM -E
	positiongodbook <- which(godbook$to == species)
		
	# here we sample until the spt is bigger than the traject length of the leaf, because the process should be valid for speciation as for extinction.
	#
	##
	while (spt <= godbook$age[positiongodbook]){
		#print(paste("again.. because spt was=", spt))
		spt <- eval(rnumbsp)
		##SM - S now we examine for shifts and update spt and shiftspm
		testshiftspout <- testshiftsp(spt)
		spt <- testshiftspout$spt
		#SM - E
	}
	## SM - S
	shiftspm <- rbind(shiftspm, testshiftspout$shiftspmr)
	#updating extinction strenght.. inheriting from previous node = "mother"
	shiftextm <- rbind(shiftextm ,c((nextsp - i),(shiftextm[shiftextm[,"node"]==species,"strength"]))) #inherit ancestor shift)
	##SM - E
	# then here we cut the overall sampled speciation time (i.e. spt) to add only the additional age
	spt <- spt - godbook$age[positiongodbook]
	##
	#
	
	
	godbook$age[positiongodbook] <-  godbook$age[positiongodbook] + spt
	#on the case of the mother, we don't sample any more extinction time, since it has already been given
	
	#here we see if there is speciation or extinction
	{
	if ( godbook$age[positiongodbook] <= godbook$destiny[positiongodbook] | is.na(godbook$destiny[positiongodbook]))
	{
		status <- "sp" #occurred an speciation
		edge.length <- c(edge.length, spt)
		leaves <- c(leaves, (nextsp - i))
		godbook$to[positiongodbook] <- nextsp -i
		ntime <- spt
	}
	else
	{
		status <- "ext" #occurred an extinction
		edge.length <- c(edge.length, godbook$destiny[positiongodbook] - sum(  trajectoryuntil(species,  godbook$from[which(godbook$to == species)])))
		extinct <- c(extinct, (nextsp - i))
		godbook$status[positiongodbook] <- status
		ntime <- edge.length[length(edge.length)]

	}
	}
	pnodges <- c(pnodges, (nextsp - i))
	timeline <- c(timeline, (ptime + ntime))
	godbook$to[positiongodbook] <- nextsp - i

	#for the daughter lineage birth
	i <- 2
	edge <- rbind( edge, c(species, (nextsp - i)))
	spt <- eval(rnumbsp)
	##SM -S now we examine for shifts and update spt and shiftspm
	testshiftspout <- testshiftsp(spt)
	spt <- testshiftspout$spt
	shiftspm <- rbind(shiftspm, testshiftspout$shiftspmr)
	##SM -E
	#to remove NaN warnings messages
	{
	if (distributionextparameters[1] == 0)
	{
		extt <- suppressWarnings(eval(rnumbext))
	}
	else
	{
		extt <- eval(rnumbext)
	}
	}
	# if to see if the user simulates with extinction = ZERO and avoid error generated by exp. distribution when rate equals zero
	{
	if (is.nan(extt))
	{
		extt <- NA # we add so that the sp will allways occurs, and we will never have an extinction (later on the NA will be replaced by the age)
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
	godbook$age <- c(godbook$age, spt)
	{
	if (spt <= extt | is.na(extt))
	{
		status <- "sp" #occurred an speciation
		edge.length <- c(edge.length, spt)
		leaves <- c(leaves, (nextsp - i))
		ntime <- spt
	}
	else
	{
		status <- "ext" #occurred an extinction
		edge.length <- c(edge.length, extt)
		extinct <- c(extinct, (nextsp - i))
		ntime <- extt
	}
	}
	godbook$status <- c(godbook$status, status)
	pnodges <- c(pnodges, (nextsp - i))
	timeline <- c(timeline, (ptime + ntime))
		
	#
	leaves <- leaves[-which(leaves == species)]
	
	#now we test for the smallest pnodge and see if there is the enough leaving creatures
	while (stopsearch == FALSE)
	{
		{
		if (length(pnodges) == 0)
		{
			stop <- TRUE	
			stopsearch <- TRUE
		}
		else
		{			
			{
			if (length(pnodges) >= m )
			{
				#we achieved m species at least
				bingo<- TRUE
				stop <- TRUE
				stopsearch <- TRUE
				age <- min(timeline)
				realleaves <- pnodges
			}
			else
			{
				positionsp <- which(timeline == min(timeline))[1]
				#testing fit the smallest tested one is an extinct
				{
				if (pnodges[positionsp] %in% extinct)
				{
					#if it is on the extinct, we will remove it from the pnodges and timeline and add it to the realextinct
					realextinct <- c(realextinct, pnodges[positionsp])
					pnodges <- pnodges[-positionsp]
					timeline <- timeline[-positionsp]
				}
				else
				{
					stopsearch <- TRUE
				}
				}
			}
			}
		}
		}	
	}
}

#after this point, all trees that are leaving this last while have at least m leaving species....
}# ######## while bingo finish

extinct <- realextinct
realleaves <- pnodges
#cutting tree at desired age
step <- 1
for (step in 1: length(pnodges))
{
	traject <- trajectory(pnodges[step])
	traject <- traject[-length(traject)]
	edge.length[which(edge[,2]==pnodges[step])] <- (age-sum(traject))
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
      	shiftedspliving[data[sort.list(data[,"node"]),"strength" ]!=1] <- 1
      	
      	tip.label[shiftedspliving==1] <- paste(tip.label[shiftedspliving==1], shiftsplabel, sep=" ")
      	shiftedspliving <- cbind(realleaves,shiftedspliving)
      	colnames(shiftedspliving) <- c("LivingSpecies", "shift")
      	#for shifts on extinction
      	shiftedextliving <- rep(0,length(realleaves))
      	
      	#shiftedextliving[shiftextm[shiftextm[,"node"]%in%realleaves,"strength"]!=1] <- 1
      	
      	data <- shiftextm[shiftextm[,"node"]%in%realleaves,]
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
	root.edge <- edge.length[1]
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
			#this is done to be able to plot and be compatible with `ape` package 
			mytree <- collapse.singles(mytree)
			#checking status of 'complete' and take or dont take extinct species out of final tree
			{
			if (complete == FALSE)	
			{
				mytree<- drop.fossil(mytree)
				mytree$root.edge <- root.edge
			}
			}		
		}
		}		
	}
	}	
}
}
return(mytree)
}
