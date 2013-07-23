\name{sim.taxa}
\alias{sim.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating General Model Trees on a fixed number of extant Taxa}
\description{
sim.taxa simulates trees on n species under a general model (i.e. any given distribution for speciation and extinction). The method allows to stop the simulation right after reaching a given number of tips/taxa or to simulate a bigger tree (m-tips) and then uniformly sample one of the n-tips trees that existed in the past until m tips. For this the gsa code from the R package TreeSim is used.
The method allows for symmetric and asymmetric binary split speciation mode. To keep all extinct tips use complete=TRUE, if complere=FALSE the tree is pruned, resulting on a so called reconstructed tree.
}
\usage{
sim.taxa(numbsim, n, m = n, distributionspname, distributionspparameters,
distributionextname = "rexp", distributionextparameters = 0, symmetric = TRUE,
complete = TRUE, labellivingsp = "sp.", labelextinctsp = "ext.", sampling = 2,
gsa = FALSE)
}
\arguments{
  \item{numbsim}{
is the number of simulated trees.
}
  \item{n}{
is the Number of tips in sampled trees (Number of extant sampled leaves).
}
  \item{m}{
is the number of standing taxa that will exist on the first generated trees, to then be sampled for n number of tips. Case gsa=FALSE, m is equal to n.
}
  \item{distributionspname}{
is the name of the desired probability function that will be used for the speciation process (e.g. distributionspname="rexp"). Note that the name should contain an `r` before it, since it refers to the randon number of the desired function (e.g. "rweibull", "runif").
}
  \item{distributionspparameters}{
are the parameters for the specific function desired for speciation. 
# IMPORTANT: this vector of function parameters must *start by the second one, since the first parameter will always be one for all the function and is added already by this function*. HINT: see the help of the desired function for more details (e.g. ?rexp) Example of parameter for a exponential distribution with lambda of one is: (distributionspparameters = c(1)). Entry in the distributionparameters can be "#", # or c(#,#) in case of more characters.
}
  \item{distributionextname}{
is the same as the distributionspname but for the probability of extinction (e.g. distributionextname = "rexp").
}
  \item{distributionextparameters}{
is the same as the distributionspparameters but for the extinction probability function. By default extinction is set to ZERO, i.e. no extinction (i.e. distributionextparameters = c(0)). Entries can be "#", # or c(#,#) in case of more characters.
}
  \item{symmetric}{
tells which macroevolutionary model should be used. If symmetric=TRUE the symmetric model will be used; if FALSE, the asymmetric model will be used. By default symmetric=TRUE.
}
  \item{complete}{
if complete = TRUE, the tree with the extinct and non-sampled lineages is returned. If complete = FALSE, the extinct and non-sampled lineages are suppressed. complete=FALSE by default.
}
  \item{labellivingsp}{
is the label that will be drawn on each tip surviving until the present. An automatic sequential number will be added to the chosen name. By default labellivingsp="sp."
}
  \item{labelextinctsp}{
is the label that will be drawn on each extinct tip. By default labelextinctsp="ext."
}
  \item{sampling}{
stochastic sampling of tips, see TreeSim for more details
}
  \item{gsa}{
gsa=TRUE indicates that the sim.gsa.taxa will be used and the n parameter will dictate the final number of tips on the tree. Note that m needs to be always bigger then n. 
If gsa = FALSE, there is no need of specifying n, once the final trees will be of size m.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\item{treearray}{Array of numbsim trees with the time since origin / most recent common ancestor being age. If tree goes extinct (only possible when mrca = FALSE), returns 0. If only one extant and no extinct tips, returns 1.}
}
\references{

O. Hagen, K. Hartmann, M. Steel and T. Stadler (2014). Age-dependent speciation explains empirical tree shape distribution. Manuscript.

T. Stadler: Simulating trees on a fixed number of extant species. Syst. Biol. (2011) 60: 676-684.
}
\author{Oskar Hagen, Tanja Stadler}
%\note{%todo
%}
\seealso{sim.age, sim.gsa.taxa}
\examples{
## example 1 ##
# without gsa, extinction and using weibull probability for speciation
my3s <- sim.taxa(1, n=10, distributionspname="rweibull", distributionspparameters=c(1.5,1),
symmetric=TRUE, complete=TRUE, labellivingsp="tip")
plot(my3s[[1]])

## example 2 ##
# using gsa=TRUE with extinction and exponential distribution for speciation and extinction
# note that this function may take a long timce since gsa=TRUE
mytree <-  sim.taxa(numbsim=10, n=10, m=15,  distributionspname="rexp",
distributionspparameters=c(1.5), distributionextname="rexp", 
distributionextparameters="0.1",symmetric = TRUE, complete=TRUE, 
labellivingsp="sp.", labelextinctsp="ext.", sampling=2, gsa=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ general model }
\keyword{ phylogenetic trees }
\keyword{ phylogenies }
\keyword{ simulation }
\keyword{ macroevolution }
