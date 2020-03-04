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

O. Hagen and T. Stadler (2017). TreeSimGM: Simulating phylogenetic trees under general Bellman Harris models with lineage-specific shifts of speciation and extinction in R. Methods in Ecology and Evolution.

O. Hagen, K. Hartmann, M. Steel and T. Stadler (2015). Age-dependent Speciation can Explain the Shape of Empirical Trees. Systematic Biology. 64, v.3, p.432-440.
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

shif_spext_asym <- sim.taxa(numbsim=1, n=5,  waitsp="rexp(1)", 
waitext="rexp(0.5)", symmetric = FALSE, 
shiftsp=list(prob=0.5, strength="rnorm(0.5,0.01)"), 
shiftext=list(prob=0.3, strength="runif(0,0.5)"))

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