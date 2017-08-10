\name{sim.age}
\alias{sim.age}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating General Model Trees on Age}
\description{
sim.age simulates stochastic trees conditioned on the time since origin (i.e. the stem age) of the extant tips under a general model (i.e. any given distribution or function for speciation and extinction). For simulations stopped when a fixed number of extant species is reached see sim.taxa.
The method allows for symmetric and asymmetric speciation mode. To keep all extinct tips in the tree use complete=TRUE, if complete=FALSE all extinct lineages without extant descendants are pruned, resulting in a so called reconstructed tree. The user needs to specify a distribution or function for the waiting time to speciation and extinction.
To account for changes in the speciation and extinction process which are inherited by the descendants, the change probability needs to be non-zero, "shiftsp$prob" and  "shiftext$prob" respectively. Since it is the probability of changes in new species of age 0, this value has to range between 0 and 1. It is possible to have changes only for the speciation or the extinction process as well as for both at the same time. By default, both speciation and extinction change probabilities are set to zero, meaning that if not stated otherwise, no changes are considered. If a change happens, a scaling factor referred to "$strength", will be multiplied to the waiting time drawn from the speciation / extinction waiting time distribution. Therefore, values smaller than 1 will shorten and values bigger than one will prolong waiting times (no negative values should be used). For every node, speciation and extinction changes are stored at "shiftsp" and "shiftext" respectively. Extant and extinct tips with changes on speciation are marked with 1 and 0 for changes or no-changes under "shifted.sp.living" and "shifted.sp.extinct" respectively. "shifted.ext.living" and "shifted.ext.extinct" follow the same order but for changes on the extinction process. Note that shifts are not allowded at the root of the tree.}
\usage{
sim.age(age, numbsim, waitsp, waitext="rexp(0)", symmetric = TRUE, complete = TRUE, 
tiplabel=c("sp.", "ext.", "Ss", "Se"), shiftsp=list(prob=0, strength="runif(0.5,0.9)"), 
shiftext=list(prob=0, strength="runif(0.1,0.2)"), sampling=list(frac=1, branchprop=FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age}{
  is the total age until each tree will be simulated (e.g. age=3), starting from one ancestral species.
  }
  \item{numbsim}{
  is the number of simulated trees.
  }
  \item{waitsp}{
  refers to the waiting time until speciation and can be informed as a (1) function that generates a single waiting time number per call or a (2) string containing the name of the probability function, i.e. random generation function and its parameters except for the first: e.g. "rexp(1.5)" or "rweibull(0.4,3)". For the last input method (2), any probability density function available in R can be used. The random number generator function has always as first parameter n (which is the number of observations). As we need only one randomly generated waiting time per draw, this first parameter of the random number generator function should be omitted from the input string. HINT: see the help of the specified distribution function for more details on its specific parameters (e.g. ?rexp). 
  }
  \item{waitext}{
  is the same as the "waitsp" but for the probability of extinction. By default, extinction is set to ZERO (e.g. waitext="rexp(0)"), i.e. no extinction process is accounted for.
  }
  \item{symmetric}{s
  defines which speciation mode should be used. If symmetric=TRUE the symmetric mode will be used; if FALSE, the asymmetric model will be used. By default symmetric=TRUE.
  }
  \item{complete}{
  if complete=TRUE, the tree with the extinct and non-sampled lineages is returned. If complete=FALSE, the extinct and non-sampled lineages are suppressed. complete=TRUE is default.
  }
  \item{tiplabel}{
  is a vector of 4 strings/labels that will be given for each tip that [1] survived until the present; [2] got extinct [3] experienced speciation changes or [4] experienced extinction changes along the way to the root. An automatic sequential number is added to the chosen label. tiplabel=c("sp.", "ext.", "Ss", "Se") is default.
  }
    \item{shiftsp}{
  a list containing [[1]] the probability $prob by which a change in the waiting time to speciation happens and the [[2]] distribution or function $strength from which a scaling factor is drawn (as in waitsp), and multiplied to the drawn speciation waiting time. shiftsp$prob should range from 0 (no change) to 1 (all new species are changed). shiftsp=list(prob=0, strength="runif(0.5,0.9)") is default.
  }
  \item{shiftext}{
  similar to shiftsp but for the shifts probabilities and strength of scaling factor of the changes in the waiting time to extinction. shiftext=list(prob=0, strength="runif(0.1,0.2)") is default
  }
  \item{sampling}{
  a list containing [[1]] the sampling fraction $frac and [[2]] a boolean $branchprop defining if the sampling should be proportional to the pendant branch lengths, in this case, longer branches would be more likely to be sampled. If $frac is smaller than 1, one sampled tree is returned with n tips as specified by the user, and the full tree prior to sampling is returned inside the tree object as $beforesampling. In the sampled tree, shift information can only be visualized though the tip labels, a complete shift history can be retrieved for the full tree $beforesampling. sampling=list(frac=1, branchprop=FALSE) is default.
  }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\item{treearray}{Array of numbsim trees with a fixed time since origin. If tree goes extinct, 0 is returned. If only one extant and no extinct tips are present, 1 is returned. For every node, including the root, speciation and extinction changes applied (i.e. the scaling factor) are stored at $shiftsp and $shiftext respectively. Extant and extinct tips with changes on speciation are marked with 1 and 0 for change or no-change under $shifted.sp.living and $shifted.sp.extinct respectively. The $shifted.ext.living and $shifted.ext.extinct follow the same logic but store shifts in the extinction process. In the case of incomplete extant species sampling, the complete simulated tree is returned via $beforesampling}
}
\references{

O. Hagen and T. Stadler. TreeSimGM: simulating phylogenetic trees under general Bellman Harris and lineage shift speciation-extinction models in R. (submitted manuscript)

O. Hagen, K. Hartmann, M. Steel and T. Stadler (2015). Age-dependent Speciation can Explain the Shape of Empirical Trees. Systematic Biology. 64, v.3, p.432-440.

T. Stadler: Simulating trees on a fixed number of extant species. Syst. Biol. (2011) 60: 676-684.
}
\author{Oskar Hagen, Tanja Stadler}
%\note{%todo
%}
\seealso{sim.taxa, track.shift}
\examples{

## example 1 ##

# Simulation of one tree with age 3 with a symmetric speciation mode, 
# a Weibull disribution for the waiting time until speciation, 
# and without extinction 

sim.age(3,1,"rweibull(0.4,3)")


## example 2 ##

# Plot one of ten simulated trees under symmetric speciation mode, with
# a uniform distribution for the waiting time to speciation and 
# an exponential distribution for the waiting time to extinction

my3s <- sim.age(age=3,numbsim=10,"runif(0.5,2)", "rexp(0.1)")

# note that for runif, the first argument is omitted and 0 stands for the
# minimum parameter of this specific function and 2 stands for the maximum

plot(my3s[[1]])


## example 3 ##

# Simulating trees with changes in the waiting time to speciation and extinction.
# The probability for a change of speciation / extinction in a new species 
# is 0.1. Upon a change in speciation, the speciation time is scaled by a factor 
# drawn from a normal distribution with mean 0.5 and sd of 0.05. This implies that 
# changes shorten the time to speciation.
# Changes in extinction are normally distributed with mean of 1.5, 
# leading to extended time to extinction

shif_spext_sym <- sim.age(age=2, numbsim=5,  waitsp="rexp(1)", waitext="rexp(0.5)", 
symmetric = TRUE, shiftsp=list(prob=0.1, strength="rnorm(0.5,0.05)"),
shiftext=list(prob=0.1, strength="rnorm(1.5,0.05)"))


## example 4 ##

# Simulating trees under an asymmetric speciation mode with changes on extinction waiting times

set.seed(10)
shif_ext_asym <- sim.age(age=2, numbsim=3, waitsp="rexp(0.8)", waitext="rexp(0.5)", 
symmetric = FALSE, shiftsp=list(prob=0.1, strength="rnorm(0.7,0.9)"), 
shiftext=list(prob=0.05, strength="runif(0.4,0.5)"))


## example 5 ##

# Simulating trees using own functions instead of strings as waiting times input
# first simulating a similar tree to example 4

set.seed(10)
shif_ext_asym_funk <- sim.age(age=2, numbsim=3, waitsp=function()rexp(1,0.8), 
waitext=function()rexp(1,0.5), 
symmetric = FALSE, shiftsp=list(prob=0.1, strength=function()rnorm(1,0.7,0.9)), 
shiftext=list(prob=0.05, strength=function()runif(1,0.4,0.5)))

# test if results are identical...
identical(shif_ext_asym, shif_ext_asym_funk)

# Now we will generate a tree based on own defined waiting time rules! only one shift 
# strength for speciation.
# In our waiting time function, we choose have exponentialy distributed waiting times 
# until speciation that are limit to be at least 0.5! 
# if they are smaller than 0.5, the waiting time will be 0.5
# remember that this function need to generate one single number.
waitspfunk <- function() {
  wt=rexp(1,1.5) 
  if(wt<0.5){wt=0.5} 
  return(wt)
}
set.seed(86)
# now we plug in our function or define it directly...
funk_tree <- sim.age(age=4, numbsim=3, waitsp=waitspfunk, 
waitext=function()rexp(1,0.9), shiftsp=list(prob=0.1, strength=function()0.5))
plot(funk_tree[[1]])


## example 6 ##

# Validation of sim.age using TreePar: 
# estimating parameters for the speciation and extinction distribution
# based on simulated trees using exponential waiting times to speciation and gamma distributed 
# (i.e. age-dependent) waiting times to extinction under an asymmetric speciation mode
\dontrun{ 
sp_la <- 3
ext_shape <- 3
ext_scale <- 2
treesTreeSimGM <- sim.age(2, 10, waitsp=paste0("rexp(",sp_la,")"), 
waitext=paste0("rgamma(",ext_shape,",",ext_scale,")"), complete=TRUE )
sptimes <- lapply(treesTreeSimGM, function(x) if  (class(x)=="phylo") getx(x))
require(TreePar) #please read TreePar documentation for installation instructions. 
** This requires Matlab or at 
** least Matlab runtime installation.
setwd("C:/YourPathToTreePar/TreePar-Matlab")
yourpath to matlab runtime
math_run_path <- 'C:/Program Files/MATLAB/MATLAB Runtime/v91/runtime/win64/' 
out <- create.mat(sptimes[[1]],path=math_run_path)
bd_out <- bd.age.optim.matlab(sptimes[[1]],path=math_run_path, sampling=1,
lambdainit=0.5,kinit=3,thetainit=0.7,
numgridpts=500)
Lcond <- "C"
lambdainit = 1
kinit = 1
thetainit = 0.5
sampling = 1
model = "G"
precision = 4
matfilename = "setup"
Param <- paste("'", lambdainit, kinit, thetainit, "'")

runCmd <- paste("sh ./run_MaxLFcn.sh ", math_run_path, matfilename, "outputML", 
               as.character(precision), Lcond, model, Param, as.character(sampling), 
               sep = " ")
system(runCmd)
}

## example 7 ##

# Validation of sim.age using TreeSim:
# trees under exponentially distributed waiting times to speciation and extinction are simulated

la=1
mu=0.5

library(TreeSim)
treesTreeSim <- sim.bd.age(2, numbsim=1000, lambda=la, mu=mu, mrca = FALSE, complete = TRUE, K = 0)
library(TreeSimGM)
treesTreeSimGM <- sim.age(2, numbsim=1000, waitsp=paste0("rexp(",la,")"),
waitext=paste0("rexp(",mu,")") )
# treesTreeSim and treesTreeSimGM have the same underlying assunptions

## compare for number of tips
# get number of tips
tipsSimTree <- unlist(lapply(treesTreeSim, function(x) if  (class(x)=="phylo") length(x[[2]])))
tipsSimTreeGM <- unlist(lapply(treesTreeSimGM, function(x) if  (class(x)=="phylo") length(x[[2]])))
# make final list
finallist <- list(SimTree=tipsSimTree, SimTreeGM=tipsSimTreeGM)
# plot
boxplot(finallist, ylab="# of tips")

## compare for oldest branching events
# get oldest branching events
branchSimTree <- unlist(lapply(treesTreeSim, function(x) if  (class(x)=="phylo") max(x[[3]])))
branchSimTreeGM <- unlist(lapply(treesTreeSimGM, function(x) if  (class(x)=="phylo") max(x[[3]])))
# make final list
finallist <- list(SimTree=branchSimTree, SimTreeGM=branchSimTreeGM)
# plot
boxplot(finallist, ylab="oldest branching events")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ general model }
\keyword{ phylogenetic trees }
\keyword{ phylogenies }
\keyword{ simulation }
\keyword{ macroevolution }
\keyword{ shift model }