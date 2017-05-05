\name{sim.taxa}
\alias{sim.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating General Model Trees on a fixed number of extant Taxa}
\description{
sim.taxa simulates trees on n extant species under a general model (i.e. any given distribution for speciation and extinction). The method allows the simulation to be stopped right after reaching a given number of co-existing tips/taxa or to simulate a bigger tree (m-tips) and then uniformly sample one of the trees with n extant tips that existed in the past prior to reaching m tips. For the later, the gsa code from the R package TreeSim is used. 
All other settings and options are equivalent to sim.age, please consult the manual for that function for details.
}
\usage{
sim.taxa(numbsim, n, m = n, distributionspname, distributionspparameters,
distributionextname = "rexp", distributionextparameters = 0, symmetric = TRUE,
complete = TRUE, labellivingsp = "sp.", labelextinctsp = "ext.", sampling = 2,
gsa=FALSE, shiftspprob=0, shiftdistributionspname="runif", 
shiftdistributionspparameters=c(0.5,0.9), shiftextprob=0, 
shiftdistributionextname="runif", shiftdistributionextparameters=c(0.1,0.2), 
shiftsplabel="Ss", shiftextlabel="Se")
}
\arguments{
  \item{numbsim}{
is the number of simulated trees.
}
  \item{n}{
is the number of tips in sampled trees (number of extant sampled leaves).
}
  \item{m}{
is the number of standing taxa that will exist on the first generated trees, to then be sampled for n number of tips. In case gsa=FALSE, m is set equal to n.

}
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
is the label that will be drawn for each tip surviving until the present. An automatic sequential number will be added to the chosen name. By default labellivingsp="sp.". For example: labellivingsp="firefly" for a simulation of a fireflies phylogeny.
}
  \item{labelextinctsp}{
is the label that will be drawn in each extinct tip. By default labelextinctsp="ext.".
}
  \item{sampling}{
Parameter determining how close the returned trees in treearray are to the "true" distribution. The higher 'sampling', the closer the output trees to the 'true' distribution.Higher values of sampling return fewer output trees meaning a larger input treearray is needed. See TreeSim::sim.gsa.taxa for more details.
}
  \item{gsa}{
gsa=TRUE indicates that the sim.gsa.taxa will be used and that the n parameter will dictate the final number of tips on the tree. Note that m needs to be always bigger then n. 
If gsa = FALSE, there is no need of specifying m, once the final trees will be of size n.
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
\seealso{sim.age, sim.gsa.taxa, track.shift}
\examples{
## example 1 ##

# Simulating trees under a Weibull distributed time to speciation, no extinction,
# and using a symmetric speciation mode. The simulation is stopped once 10 co-existing species
# exist for the first time (gsa = FALSE).

my3 <- sim.taxa(1, n=10, distributionspname="rweibull", distributionspparameters=c(1.5,1),
symmetric=TRUE, complete=TRUE, labellivingsp="tip")

plot(my3[[1]])


## example 2 ##

# Trees are simulated using gsa=TRUE with an exponential distribution for 
# speciation and extinction

\dontrun{ 
** long runing timce since gsa=TRUE **
mytree <-  sim.taxa(numbsim=10, n=10, m=15,  distributionspname="rexp",
distributionspparameters=c(1.5), distributionextname="rexp", 
distributionextparameters=c(0.1),symmetric = TRUE, complete=TRUE, 
labellivingsp="sp.", labelextinctsp="ext.", sampling=2, gsa=TRUE)
}

## example 3 ##

# Now changes in speciation waiting times occur in new species with 
# probability 0.1, again under a symmetric mode

shift_sp_sym <- sim.taxa(numbsim=10, n=10, distributionspname="rexp", 
distributionspparameters=c(1), distributionextname="rexp", 
distributionextparameters=c(0.5),symmetric = TRUE, shiftspprob=0.1, 
shiftdistributionspname="runif",  shiftdistributionspparameters=c(0.5,0.9))


## example 4 ##

# Simulations under an asymmetric speciation mode with changes in 
# speciation and extinction waiting times

shif_spext_asym <- sim.taxa(numbsim=10, n=10, distributionspname="rexp", 
distributionspparameters=c(1), distributionextname="rexp", 
distributionextparameters=c(0.5),symmetric = FALSE, shiftspprob=0.1, 
shiftdistributionspname="runif", shiftdistributionspparameters=c(0.5,0.9), 
shiftextprob=0.05, shiftdistributionextname="runif", 
shiftdistributionextparameters=c(0,0.5))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ general model }
\keyword{ phylogenetic trees }
\keyword{ phylogenies }
\keyword{ simulation }
\keyword{ macroevolution }
\keyword{ shift model }