\name{sim.age}
\alias{sim.age}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating General Model Trees on Age}
\description{
sim.age simulates stochastic trees conditioned on the time since origin (i.e. the stem age) of the extant tips under a general model (i.e. any given distribution for speciation and extinction). For simulations stopped once a fixed number of extant species is reached see sim.taxa.
The method allows for symmetric and asymmetric speciation mode. To keep all extinct tips in the tree use complete=TRUE, if complete=FALSE all extinct lineages without extant descendants are pruned, resulting in a so called reconstructed tree. The user needs to specify a distribution for the waiting time to speciation and extinction.
To account for changes in the speciation and extinction process which are inherited by the descendants, the change probability needs to be non-zero, "shiftspprob" and  "shiftextprob" respectively. Since it is the probability of changes in new species of age 0, this value has to range between 0 and 1. It is possible to have changes only for the speciation or the extinction process as well as for both at the same time. By default, both speciation and extinction change probabilities are set to zero, meaning that if not stated otherwise, no changes are considered. If a change happens, a scaling factor referred to "shiftstrength", will be multiplied to the waiting time drawn from the speciation / extinction waiting time distribution. Therefore, values smaller than 1 will shorten and values bigger than one will prolong waiting times (no negative values should be used). For every node, speciation and extinction changes are stored at "shiftsp" and "shiftext" respectively. Extant and extinct tips with changes on speciation are marked with 1 and 0 for changes or no-changes under "shifted.sp.living" and "shifted.sp.extinct" respectively. "shifted.ext.living" and "shifted.ext.extinct" follow the same order but for changes on the extinction process. Note that shifts are not allowded at the root of the tree.}
\usage{
sim.age(age, numbsim, distributionspname, distributionspparameters,
distributionextname = "rexp", distributionextparameters = 0, symmetric = TRUE, 
complete = TRUE, labellivingsp = "sp.", labelextinctsp = "ext.", shiftspprob=0, 
shiftdistributionspname="runif", shiftdistributionspparameters=c(0.5,0.9), 
shiftextprob=0, shiftdistributionextname="runif", 
shiftdistributionextparameters=c(0.1,0.2), shiftsplabel="Ss", shiftextlabel="Se")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age}{is the total age until each tree will be simulated (e.g. age=3), starting from one ancestral species.}
  \item{numbsim}{is the number of simulated trees.}
  \item{distributionspname}{
is the name of the probability function that will be used for the speciation process (e.g. distributionspname="rexp"). Any probability density function available in R can be used. Note that by R convention, the name has to contain an r before it, since it refers to the randon number of the desired function (e.g. "rweibull", "runif").
  }
  \item{distributionspparameters}{
are the parameters for the specific speciation probability density function. 
# IMPORTANT: the first input of the random number generator function specified in distributionspname (e.g. "rexp") is the number of draws from that distribution. As we always draw one value, this entry needs to be ignored, and the input for distribuionspparameters starts with the second input of say "rexp". HINT: see the help of the specified distribuitonspparameters function for more details (e.g. ?rexp). Example of parameter for a exponential distribution with parameter of one is: (distributionspparameters = c(1)).
}
  \item{distributionextname}{
is the same as the "distributionspname" but for the probability of extinction (e.g. distributionextname = "rexp").
}
  \item{distributionextparameters}{
is the same as the "distributionspparameters" but for the extinction probability function. By default extinction is set to ZERO, i.e. no extinction.
}
  \item{symmetric}{
defines which speciation mode should be used. If symmetric=TRUE the symmetric mode will be used; if FALSE, the asymmetric model will be used. By default symmetric=TRUE.
}
  \item{complete}{
if complete=TRUE, the tree with the extinct and non-sampled lineages is returned. If complete=FALSE, the extinct and non-sampled lineages are suppressed. complete=TRUE is default.
}
  \item{labellivingsp}{
is the label that will be drawn for each tip surviving until the present. An automatic sequential number will be added to the chosen name. By default labellivingsp="sp.". For example: labellivingsp="bird" for a simulation of a birds phylogeny.
}
  \item{labelextinctsp}{
is the label that will be drawn in each extinct tip. By default labelextinctsp="ext.".
}
    \item{shiftspprob}{
is the probability by which a change in the waiting time to speciation happens, default is 0, that means no change. This value should range from 0 (no change) to 1 (all new species are changed).
}
  \item{shiftdistributionspname}{
if a speciation shift happens, "shiftdistributionspname" is the distribution from which a scaling factor is drawn, which is multiplied to the drawn speciation waiting time.
}
  \item{shiftdistributionspparameters}{
are the parameters of the chosen distribution in "shiftdistributionspname". To be defined as parameters in distributionspparameters.
}
  \item{shiftextprob}{
is the probability by which a change in waiting time to extinction happens, default is 0, that means no change. This value should range from 0 (no change) to 1 (all species are changed).
}
  \item{shiftdistributionextname}{
if a extinction shift happens, shiftdistributionspname is the distribution from which a scaling factor is drawn, which is multiplied to the drawn speciation waiting time.
}
  \item{shiftdistributionextparameters}{
are the parameters of the chosen distribution in "shiftdistributionextname". To be defined as parameters in distributionextparameters.
}
  \item{shiftsplabel}{
label to be added to the tips that experienced speciation changes along the way to the root.
}
  \item{shiftextlabel}{
label to be added to the tips that experienced extinction changes along the way to the root.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\item{treearray}{Array of numbsim trees with a fixed time since origin. If tree goes extinct, 0 is returned. If only one extant and no extinct tips are present, 1 is returned. For every node, including the root, speciation and extinction changes (i.e. the scaling factor) are stored at $shiftsp and $shiftext respectively. Extant and extinct tips with changes on speciation are marked with 1 and 0 for change or no-change under $shifted.sp.living and $shifted.sp.extinct respectively. The $shifted.ext.living and $shifted.ext.extinct follow the same order but store shifts in the extinction process.}
}
\references{

O. Hagen and T. Stadler. TreeSimGM: simulating phylogenetic trees under general speciation-extinction models in R. (submitted manuscript)

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

sim.age(3,1,"rweibull", c(0.4,3))


## example 2 ##

# Plot one of ten simulated trees under symmetric speciation mode, with
# a uniform distribution for the waiting time to speciation and 
# an exponential distribution for the waiting time to extinction

my3s <- sim.age(age=3,numbsim=10,"runif", c(0,2), "rexp", c(0.1) )

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

shif_spext_sym <- sim.age(age=2, numbsim=5,  distributionspname="rexp", 
distributionspparameters=c(1), distributionextname="rexp", 
distributionextparameters=c(0.5),symmetric = TRUE, shiftspprob=0.1, 
shiftdistributionspname="rnorm", shiftdistributionspparameters=c(0.5,0.05),
shiftextprob=0.1, shiftdistributionextname="rnorm", 
shiftdistributionextparameters=c(1.5,0.05))


## example 4 ##

# Simulating trees under an asymmetric speciation mode with changes on extinction waiting times

set.seed(10)
shif_ext_asym <- sim.age(age=2, numbsim=5,distributionspname="rexp", 
distributionspparameters=c(1), distributionextname="rexp", distributionextparameters=c(0.5),
symmetric = FALSE, labelextinctsp="ext.", shiftspprob=0.1, shiftdistributionspname="runif", 
shiftdistributionspparameters=c(0.5,0.9), shiftextprob=0.05, shiftdistributionextname="runif", 
shiftdistributionextparameters=c(0.3,0.5))


## example 5 ##

# Validation of sim.age using TreePar: 
# estimating parameters for the speciation and extinction distribution
# based on simulated trees using exponential waiting times to speciation and gamma distributed 
# (i.e. age-dependent) waiting times to extinction under an asymmetric speciation mode
\dontrun{ 
sp_la <- 3
ext_shape <- 3
ext_scale <- 2
treesTreeSimGM <- sim.age(2, 10, distributionspname="rexp", 
distributionspparameters=c(sp_la), distributionextname="rgamma", 
distributionextparameters=c(ext_shape,ext_scale), complete=TRUE )
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

## example 6 ##

# Validation of sim.age using TreeSim:
# trees under exponentially distributed waiting times to speciation and extinction are simulated

la=1
mu=0.5

library(TreeSim)
treesTreeSim <- sim.bd.age(2, numbsim=1000, lambda=la, mu=mu, mrca = FALSE, complete = TRUE, K = 0)
library(TreeSimGM)
treesTreeSimGM <- sim.age(2, numbsim=1000, distributionspname="rexp", distributionspparameters=la,
distributionextname="rexp", distributionextparameters=mu )
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