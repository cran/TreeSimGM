\name{track.shift}
\alias{track.shift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tracking back shifts (speciation or extinction) on TreeSimGM simulated tree}
\description{
track.shift returns all changes in speciation or extinction between the root of a tree and a specific node or tip in a phylogenetic tree simulated with TreeSimGM. It can track back changes on the speciation process (shift='sp') or on the extinction process (shift='ext'). Changes are the multipliers of the drawn waiting times.
}
\usage{
track.shift(shift="sp",tree, node)
}
\arguments{
  \item{shift}{
specifies which macroevolutionary process shall be tracked. Can be "sp" or "ext" for tracking speciation or extinction shifts repectively (default is "sp").
}
  \item{tree}{
is the simulated tree object.
}
  \item{node}{
is the node number or tip label that should be tracked. The numbers are the number of the attributes tip.label or as specified in '$edge' on the phylo tree object.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
\item{out}{is a vector of changes, i.e. scaling factors of the waiting times process, from MRCA until the desired node. Each reported scaling factor is named with a node number (as in $edge). The scaling factor change happened at the start of the branch ancestral to that named node.}
}
\references{

O. Hagen and T. Stadler. TreeSimGM: simulating phylogenetic trees under general speciation-extinction models in R. (submitted manuscript)

}
\author{Oskar Hagen, Tanja Stadler}
%\note{%todo
%}
\seealso{sim.age, sim.taxa}
\examples{
# setting seed for reproductibility

set.seed(2142) 

# Simulating a tree under an asymmetric speciation mode with frequent changes in 
# speciation and extinction (repectively 50% and 30% probability)

shif_spext_asym <- sim.taxa(numbsim=1, n=5,  distributionspname="rexp", 
distributionspparameters=c(1), distributionextname="rexp", 
distributionextparameters=c(0.5),symmetric = FALSE, shiftspprob=0.5, 
shiftdistributionspname="rnorm", shiftdistributionspparameters=c(0.5,0.01), 
shiftextprob=0.3, shiftdistributionextname="runif", 
shiftdistributionextparameters=c(0,0.5))

# visualiying the generated stochastic tree

plot(shif_spext_asym[[1]]) 


# tracking speciation shifts for species 3
# note that this function returns the shifts on the branches
# starting from the mrca

track.shift(shift="sp", tree=shif_spext_asym[[1]], node=3)

# tracking speciation shifts for species 7

track.shift(shift="sp", tree=shif_spext_asym[[1]], node=7)

# tracking extinction shifts for species 7

track.shift(shift="ext", tree=shif_spext_asym[[1]], node=7)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ general model }
\keyword{ phylogenetic trees }
\keyword{ phylogenies }
\keyword{ simulation }
\keyword{ macroevolution }
\keyword{ shift model }