\name{sim.age}
\alias{sim.age}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating General Model Trees on Age}
\description{
sim.age simulates stochastic trees conditioned on the time since origin of the extant tips under a general model (i.e. any given distribution for speciation and extinction). The method allows for symmetric and asymmetric binary split speciation mode. To keep all extinct trees use complete=TRUE, if complete=FALSE the tree is pruned, resulting on a so called reconstructed tree.}
\usage{
sim.age(age, numbsim, distributionspname, distributionspparameters, distributionextname = "rexp", distributionextparameters = 0, symmetric = TRUE, complete = TRUE, labellivingsp = "sp.", labelextinctsp = "ext.")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age}{is the total age until each tree will be simulated (e.g. age=3). Time since origin / most recent common ancestor}
  \item{numbsim}{is the number of simulated trees}
  \item{distributionspname}{is the name of the desired probability function that will be used for the speciation process (e.g. distributionspname="rexp"). Note that the name should contain an `r` before it, since it refers to the randon number of the desired function (e.g. "rweibull", "runif")}
  \item{distributionspparameters}{are the parameters for the specific function desired for speciation. 
# IMPORTANT: this vector of function parameters must *start by the second one, since the first parameter will always be one for all the function and is added already by this function*. HINT: see the help of the desired function for more details (e.g. ?rexp) Example of parameter for a exponential distribution with lambda of one is: (distributionspparameters = c(1)). Entry in the distributionparameters can be "#", # or c(#,#) in case of more characters}
  \item{distributionextname}{is the same as the distributionspname but for the probability of extinction (e.g. distributionextname = "rexp")}
  \item{distributionextparameters}{is the same as the distributionspparameters but for the extinction probability function. By default extinction is set to ZERO, i.e. no extinction (i.e. distributionextparameters = c(0)). Entries can be "#", # or c(#,#) in case of more characters}
  \item{symmetric}{tells which macroevolutionary model should be used. If symmetric=TRUE the symmetric model will be used; if FALSE, the asymmetric model will be used. By default symmetric=TRUE}
  \item{complete}{if complete = TRUE, the tree with the extinct and non-sampled lineages is returned. If complete = FALSE, the extinct and non-sampled lineages are suppressed. complete=FALSE by default}
  \item{labellivingsp}{is the label that will be drawn on each tip surviving until the present. An automatic sequential number will be added to the chosen name. By default labellivingsp="sp."}
  \item{labelextinctsp}{is the label that will be drawn on each extinct tip. By default labelextinctsp="ext."}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\item{treearray}{Array of numbsim trees with the time since origin / most recent common ancestor being age. If tree goes extinct (only possible when mrca = FALSE), returns 0. If only one extant and no extinct tips, returns 1.}
}
\references{
O. Hagen, K. Hartmann, M. Steel and T. Stadler (2013). Age-dependent speciation explains empirical tree shape distribution. Manuscript.

T. Stadler: Simulating trees on a fixed number of extant species. Syst. Biol. (2011) 60: 676-684.
}
\author{Oskar Hagen, Tanja Stadler}
%\note{%todo
%}
\seealso{sim.taxa}
\examples{
# One tree with age 3 without extinction and weibull 
# disribution for the waiting time until speciation
sim.age(3,1,"rweibull", c(0.4,3))

# Plot one of ten trees
# uniform distribution for speciation and 
# exponential distribution for extinction
xuxu <- sim.age(age=3,10,"runif", c(0,2), "rexp", c(0.1) )
# note that for runif, the first argument is omitted and 0 stands for the
# min parameter of this specific function as 2 stands for max.
plot(xuxu[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ general model }
\keyword{ phylogenetic trees }
\keyword{ phylogenies }
\keyword{ simulation }
\keyword{ macroevolution }

