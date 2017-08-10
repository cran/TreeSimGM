sim.taxa <-
function (numbsim, n, m=n,  waitsp, waitext="rexp(0)", symmetric=TRUE, complete=TRUE, tiplabel=c("sp.", "ext.","Ss", "Se"), 
          shiftsp=list(prob=0, strength="runif(0.5,0.9)"), shiftext=list(prob=0, strength="runif(0.1,0.2)"), 
          sampling=list(frac=1, branchprop=FALSE), sampling.gsa=1, gsa=FALSE) {
# numbsim is the number of simulated trees
# n is the Number of tips in sampled trees (Number of extant sampled leaves)
# m is the number of standing taxa that will exist on the first generated trees, to then be sampled for n number of tips. Case gsa=TRUE, m is equal to n.
# distributionspname is the name of the desired probability function that will be used for the speciation process (e.g. distributionspname <- "rexp"). Note that the name should contain an `r` before it, since it refers to the random number of the desired function (e.g. "rweibull", "runif")
# distributionspparameters are the parameters for the specific function desired for speciation. 
# IMPORTANT: this vector of fuction parameters must *start by the second one, since the first parameter will always be one for all the function and is added already by this function*. HINT: see the help of the desired function for more details (e.g. ?rexp) Example of parameter for a exponential distribution with lambda of one (distributionspparameters <- c(1)). Entry in the distributionparameters can be "#", # or c(#,#) in case of more characters
# distributionextname is the same as the distributionspname but for the probability of extinction (e.g. distributionextname <- "rexp") 
# distributionextparameters is the same as the distributionspparameters but for the extinction probability function. By default extinction is set to ZERO, i.e. no extinction (e.g. distributionextparameters <- c(0)). Entry in can be "#", # or c(#,#) in case of more characters
# symmetric tells which macro-evolutionary model should be used. If symmetric=TRUE the symmetric model will be used, else if FALSE, asymmetric model will be used. By default symmetric=TRUE
# complete: If complete = TRUE, the tree with the extinct and non-sampled lineages is returned. If complete = FALSE, the extinct and non-sampled lineages are suppressed. Complete=FALSE by default
# labellivingsp is the label that will be drawn on each tip surving until the present. An automatic sequential number will be added to the chosen name. By default labellivingsp="sp."
# labelextinctsp is the label that will be drawn on each extinct tip. By default labelextinctsp <- "ext."
# sampling: stochastic sampling, default
# gsa TRUE indicates that the sim.gsa.taxa will be used, the n parameter indicates the final number of species. Note that m needs to be allways bigger then n. If gsa = FALSE, there is no need of specifying n, once the final trees will be of size m
# entry in the distributionparameters can be "#", # or c(#,#) in case of more variables
# shiftspprob: frequency by which a speciation shift happens, default is 0, that means no shift. This value should range from 0 (no shift) to 1 (all species are shifted)
# shiftdistributionspname: distribution by which the shift (waiting time multiplier) will be drawn
# shiftdistributionspparameters: parameters of the chosen distribution  
# shiftextprob: frequency by which a extinction shift happens, default is 0, that means no shift. This value should range from 0 (no shift) to 1 (all species are shifted)
# shiftdistributionextname: distribution by which the shift (waiting time multiplier) will be drawn
# shiftdistributionextparameters: parameters of the chosen distribution 
# shiftsplabel: label to be added to the species that suffered speciation shift
# shiftextlabel: label to be added to the species that suffered extinction shift
  
  if (complete==TRUE & sampling$frac!=1) {
    warning("Sampling on taxa based can only be used with complete=FALSE, thus complete was changed to FALSE")
    complete=FALSE
  }
  
  if (sampling$frac>1 | sampling$frac<0){
    warning("Sampling Sampling fraction needs to range between 0 and 1, thus sampling$frac was changed to 1")
    sampling$frac=1
  }
  
  if (sampling$frac!=1){
    n=round(n/sampling$frac)
    if (m<n) {
      warning("You are using sampling, thus tips=n/sampling$frac. m is smaller than n, thus we changed m=n/frac")
      m=n
    }
  }
  
  if (m<n){
    warning("m can not be samller than n, thus we changed m=n")
    m=n
  }
  
	check<-gsa
	if (gsa==F && complete==T){check<-T}
	mytreegsazed <- list()
  while (length(mytreegsazed) < numbsim)
  {
  	mytree <- list()
  	step <- 1
  	{
  	if (symmetric == TRUE) 	
  	{
  		for (step in 1: (numbsim) ){
  			mytreenext <- mytree.symmetric.taxa(m=m, waitsp=waitsp, waitext=waitext, complete=check, tiplabel=tiplabel, 
  			                                    shiftsp=shiftsp, shiftext=shiftext, sampling=sampling, gsa=gsa)
  			mytree<- c(mytree, list(mytreenext))
  		}	
  	}
  	else
  	{
  		for (step in 1: (numbsim) ){
  			mytreenext <- mytree.asymmetric.taxa(m=m, waitsp=waitsp, waitext=waitext, complete=check, tiplabel=tiplabel,
  			                                     shiftsp=shiftsp, shiftext=shiftext, sampling=sampling, gsa=gsa)
  			mytree<- c(mytree, list(mytreenext))
   		}
  	}
  	}
  	{
  	if (gsa==T)
  	{
  		mytreegsa <- sim.gsa.taxa(mytree, n=n, sampling=sampling.gsa, frac=sampling$frac, complete=complete)
  	} 
  	else 
  	{
  		mytreegsa <- mytree #gsa is FALSE
  	}
  	}
  	mytreegsazed <- c(mytreegsazed, mytreegsa)
	}
mytreegsazeds <- sample(mytreegsazed, numbsim)
mytreegsazed <- mytreegsazeds
return(mytreegsazed)	
}
