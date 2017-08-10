\name{sim.taxa}
\alias{sim.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating General Model Trees on a fixed number of extant Taxa}
\description{
sim.taxa simulates trees on n extant species under a general model (i.e. any given distribution for speciation and extinction). The method allows the simulation to be stopped right after reaching a given number of co-existing tips/taxa or to simulate a bigger tree (m-tips) and then uniformly sample one of the trees with n extant tips that existed in the past prior to reaching m tips. For the later, the gsa code from the R package TreeSim is used. 
All other settings and options are equivalent to sim.age, please consult the manual for that function for details.
}
\usage{
sim.taxa(numbsim, n, m = n, waitsp, waitext="rexp(0)",symmetric = TRUE, complete = TRUE, 
tiplabel=c("sp.", "ext.", "Ss", "Se"), shiftsp=list(prob=0, strength="runif(0.5,0.9)"), 
shiftext=list(prob=0, strength="runif(0.1,0.2)"), sampling=list(frac=1, branchprop=FALSE),
sampling.gsa=1, gsa=FALSE)
}
\arguments{
  \item{numbsim}{
is the number of simulated trees.
}
  \item{n}{
is the number of tips in the sampled trees (number of extant sampled leaves).
}
  \item{m}{
is the number of standing taxa that will exist on the first generated trees, to then be sampled for n number of tips. In case gsa=FALSE, m is set equal to n.
}
\item{waitsp}{
  refers to the waiting time until speciation and can be informed as a (1) function that generates a single waiting time number per call or a (2) string containing the name of the probability function, i.e. random generation function and its parameters except for the first: e.g. "rexp(1.5)" or "rweibull(0.4,3)". For the last input method (2), any probability density function available in R can be used. The random number generator function has always as first parameter n (which is the number of observations). As we need only one randomly generated waiting time per draw, this first parameter of the random number generator function should be omitted from the input string. HINT: see the help of the specified distribution function for more details on its specific parameters (e.g. ?rexp). 
}
  \item{waitext}{
  is the same as the "waitsp" but for the probability of extinction. By default, extinction is set to ZERO (e.g. waitext="rexp(0)"), i.e. no extinction process is accounted for.
}
  \item{symmetric}{
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
  a list containing [[1]] the sampling fraction $frac and [[2]] a boolean $branchprop defining if the sampling should be proportional to the pendant branch lengths, in this case, longer branches would be more likely to be sampled. If $frac is smaller than 1, one sampled tree is returned with n tips as specified by the user, and the full tree is returned inside the tree object as $beforesampling. In the sampled tree, shift information can only be visualized though the tip labels, a complete shift history can be retrieved at the full tree $beforesampling. sampling=list(frac=1, branchprop=FALSE) is default.
}
  \item{sampling.gsa}{
Parameter determining how close the returned trees in treearray are to the "true" distribution. The higher 'sampling', the closer the output trees to the 'true' distribution.Higher values of sampling return fewer output trees meaning a larger input treearray is needed. See TreeSim::sim.gsa.taxa for more details.
}
  \item{gsa}{
gsa=TRUE indicates that the sim.gsa.taxa will be used and that the n parameter will dictate the final number of tips on the tree. Note that m needs to be always bigger then n. 
If gsa = FALSE, there is no need of specifying m, once the final trees will be of size n.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\item{treearray}{array of numbsim trees with a fixed number of living tips. For every node, including the root, speciation and extinction changes applied (i.e. the scaling factor) are stored at $shiftsp and $shiftext respectively. Extant and extinct tips with changes on speciation are marked with 1 and 0 for change or no-change under $shifted.sp.living and $shifted.sp.extinct respectively. The $shifted.ext.living and $shifted.ext.extinct follow the same logic but store shifts in the extinction process. In the case of incomplete extant species sampling, the complete simulated tree is returned via $beforesampling}
}
\references{

O. Hagen and T. Stadler. TreeSimGM: simulating phylogenetic trees under general Bellman Harris and lineage shift speciation-extinction models in R. (submitted manuscript)

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

my3 <- sim.taxa(1, n=10, waitsp="rweibull(1.5,1)",
symmetric=TRUE, complete=TRUE, tiplabel=c("tip", "tip" ,"", ""))

plot(my3[[1]])


## example 2 ##

# Trees are simulated using gsa=TRUE with an exponential distribution for 
# speciation and extinction

\dontrun{ 
** long runing timce since gsa=TRUE **
mytree <-  sim.taxa(numbsim=10, n=10, m=15,  waitsp="rexp(1.5)",
waitext="rexp(0.1)",symmetric = TRUE, complete=TRUE, 
tiplabel=c("sp.", "ext.", NA, NA), sampling.gsa=2, gsa=TRUE)
}

## example 3 ##

# Now changes in speciation waiting times occur in new species with 
# probability 0.1, again under a symmetric mode

shift_sp_sym <- sim.taxa(numbsim=10, n=10, waitsp="rexp(1)", 
waitext="rexp(0.5)", symmetric = TRUE, shiftsp=list(prob=0.1, strength="runif(0.5,0.9)"))


## example 4 ##

# Simulations under an asymmetric speciation mode with changes in 
# speciation and extinction waiting times

shif_spext_asym <- sim.taxa(numbsim=10, n=10, waitsp="rexp(1)", 
waitext="rexp(0.5)", symmetric = FALSE, shiftsp=list(prob=0.1, strength="runif(0.5,0.9)"), 
shiftext=list(prob=0.05, strength="runif(0,0.5)"))

## example 5 ##

# Waiting times as functions instead of strings, allowing for more flexibility
t1 <- sim.taxa(1,8, waitsp="rnorm(0.5,0)")
plot(t1[[1]])
t2 <- sim.taxa(1,8, waitsp=function() rnorm(1,0.5,0))
plot(t2[[1]])
identical(t1,t2)
t3 <- sim.taxa(1,8, waitsp=function() 0.5)
plot(t3[[1]])
identical(t1,t3)
# all implementations are identical!

# Creating own function for waitsp and shiftsp
waitspfunk <- function() {
  wt=rexp(1,1.5) 
  if(wt>1.5){wt=20} 
  return(wt)
}
# here we force all values bigger than 1.5 to be very large, i.e. 20
set.seed(13)
tshittfunk <- sim.taxa(1,10, waitsp=waitspfunk, 
shiftsp=list(prob=0.2, strength=function()sample(c(0.05,0.5,0.9), 1) ))
plot(tshittfunk [[1]])



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ general model }
\keyword{ phylogenetic trees }
\keyword{ phylogenies }
\keyword{ simulation }
\keyword{ macroevolution }
\keyword{ shift model }