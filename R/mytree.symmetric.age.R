mytree.symmetric.age <-
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
  

  
# # # # # # # # # DECLARATIONS  MACHINE # # # # # # # # # 
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


  	
  #testing for shift speciation
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
  
  
  #testing for shifts on extinction process
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
  ## SM - E
  
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

##SM - S creating shiftsp and shiftext matrixes...
shiftspm <- matrix(c(-1,1,-2,1), byrow=TRUE, ncol=2, dimnames=list(NULL,c("node","strength")))
shiftextm <- matrix(c(-1,1,-2,1), byrow=TRUE, ncol=2, dimnames=list(NULL,c("node","strength")))
##SM - E

# initial if, in case the (-1,-2) edge get extinct or bigger than age
spt <- rnumbsp()#eval(rnumbsp)

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

# to see if the user simulates with extinction = ZERO and avoid error generated by expression-distribution when rate equals zero
{
if (is.nan(extt))
{
	extt <- spt +1 # we add so that the sp will allways occurs first, and we will never have an extinction
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
# for the while (increase of the tree)
while (stop == FALSE) 
{
	species <- leaves[1]
	nextsp <- min(edge[,2])
	i <- 1
		for (i in 1:2)
		{ 
			edge <- rbind( edge, c(species, (nextsp - i) ) )
			
			spt <- rnumbsp()#eval(rnumbsp)
			
			##SM - S now we examine for shifts and update spt and shiftspm
			testshiftspout <- testshiftsp(spt)
			spt <- testshiftspout$spt
			shiftspm <- rbind(shiftspm, testshiftspout$shiftspmr)
			##SM - E
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
				extt <- spt +10000 # we add this so that the sp will allways occurs first, so that we will never have an extinction
			}
			}
			##SM - S now we examine for shifts and update extt and shiftextm
			testshiftextout <- testshiftext(extt)
			extt <- testshiftextout$extt
			shiftextm <- rbind(shiftextm, testshiftextout$shiftextmr)
			##SM - E
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
					realleaves <- c( realleaves, leaves[length(leaves)] )
					leaves <- leaves[-length(leaves)]					
				}
				else
				{
					realleaves <- c( realleaves, extinct[length(extinct)] )			
					extinct <- extinct[-length(extinct)]
				}
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
						## SM - S changing nodes names
						shiftspm[ which(shiftspm[,"node"] == prealleaves[i]), "node" ] <- realleaves[i]
						shiftextm[ which(shiftextm[,"node"] == prealleaves[i]), "node" ] <- realleaves[i]
						## SM - E
				}
      	tip.label <- paste(labellivingsp, realleaves, sep = "")
      	## SM - S creating vector with shifted species and updating shifted names
      	#for shifts on speciation
      	shiftedspliving <- rep(0,length(realleaves)) #preparing vector
      	data <- shiftspm[shiftspm[,"node"]%in%realleaves,]
      	#error handling in case only one species got extinct and data is transformed into a vector and not a matrix...
      	if(class(data)=="numeric"){
      	  data <- matrix(data, ncol=2)
      	}
      	shiftedspliving[data[sort.list(data[,1]),2]!=1] <- 1
      	tip.label[shiftedspliving==1] <- paste(tip.label[shiftedspliving==1], shiftsplabel, sep=" ")
      	shiftedspliving <- cbind(realleaves,shiftedspliving)
      	colnames(shiftedspliving) <- c("LivingSpecies", "shift")
      	
      	#for shifts on extinction
      	shiftedextliving <- rep(0,length(realleaves))
      	data <- shiftextm[shiftextm[,"node"]%in%realleaves,]
      	#error handling in case only one species got extinct and data is transformed into a vector and not a matrix...
      	if(class(data)=="numeric"){
      	  data <- matrix(data, ncol=2)
      	}
      	shiftedextliving[data[sort.list(data[,1]),2]!=1] <- 1
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
  			## SM - S changing nodes names
  			shiftspm[ which(shiftspm[,"node"] == pextinct[i]), "node" ] <- extinct[i]
  			shiftextm[ which(shiftextm[,"node"] == pextinct[i]), "node" ] <- extinct[i]
  			## SM - E
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
		#substituting...
		i <- 1
		for (i in 1:length(potheredges))
		{
			edge[ edge == potheredges[i] ] <- otheredges[i]
			## SM - S changing nodes names
			shiftspm[shiftspm[,"node"]==potheredges[i],"node"] <- otheredges[i]
			shiftextm[shiftextm[,"node"]==potheredges[i],"node"] <- otheredges[i]
			## SM - E
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
		#in case only one specie is surviving, even if other specification events occurred in the history
		mytree <- 1
	}
	else
	{
		{
		if (length(realleaves)==1 & length(extinct)==0 & complete==TRUE)
		{
			#in case only one species is surviving and was the only existing one
			mytree <- 1
		}
		else
		{
			#in case non of the above condition is fulfilled, there will be a tree with no initial branch, tree starts at the MRCA
			#this is done to be able to plot using the `ape` package 
			mytree <- collapse.singles(mytree)
			#cheking status of 'complete'. Take or dont take extincted species out of final tree
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
