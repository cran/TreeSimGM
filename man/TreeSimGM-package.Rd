\name{TreeSimGM-package}
\alias{TreeSimGM-package}
\alias{TreeSimGM}
\docType{package}
\title{Simulating Phylogenies under a General Model (symmetric & asymmetric)}
\description{
TreeSimGM is a flexible simulation tool for phylogenetic trees that uses a general model. This model can assume any probability distribution for the waiting time until speciation, as well as any distribution for the waiting time until extinction. In clear English this means that 'TreeSimGM, allows you to simulate stochastic phylogenetic trees using any probability distribution and parameters for speciation and  extinction.Moreover, a general model can be distinguished in the way that speciation and extinction processes happen. In this package, we provide two general binary splitting models: (i) symmetric, in which for every speciation event new waiting times until speciation and extinction are drawn for both daughter lineages; and (ii) asymmetric, in which a speciation event results in one species with new waiting times, and another that carries the extinction time and age of its ancestor. Those models were inspired by allopatric and peripatric speciation modes respectively. Both models (symmetric and asymmetric) were created and implemented so that different processes (distributions) for speciation and extinction could be independently and explicitly specified. It is also possible to have an implicit extinction process by setting the extinction rate to zero and simulating only based on speciation processes.
}
\details{
\tabular{ll}{
Package: \tab TreeSimGM\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-02-20\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
}
TreeSimGM can simulate many trees returning a list of phylogenetic trees. The trees can be simulated in different ways: (i) by age `sim.age` (i.e. the user specifies a time, or tree age, when the simulation should stop); (ii) by taxa `sim.taxa` (i.e. the user specifies the number of tips desired on the final tree) Moreover, the user needs to specify the specialization mode (i.e. symmetric or asymmetric), thus, there are four major possible combinations: (a) sim.age with parameter symmetric = TRUE; (b) sim.age with parameter symmetric = FALSE (i.e. asymmetric); (c) sim.taxa with parameter symmetric = TRUE; (d) sim.taxa with parameter symmetric = FALSE (i.e. asymmetric)
}
\author{
Oskar Hagen, Tanja Stadler

Maintainer: Oskar Hagen <oskar.is.hagen@gmail.com>
}
\references{
#to cite this package use:  citation("TreeSimGM")

O. Hagen, K. Hartmann, M. Steel and T. Stadler (2013). Age-dependent speciation explains empirical tree shape distribution. Manuscript.

T. Stadler: Simulating trees on a fixed number of extant species. Syst. Biol. (2011) 60: 676-684.
}
%\keyword{tree, phylogeny, simulation, general model, macroevolution }
\seealso{
\code{\link[TreeSim]{TreeSim}}
\code{\link[ape]{ape}}
}
\examples{
##plots the first tree of a list of three generated
#simulation based on age. Note that by default, symmetric = TRUE
trialtheetrees <- sim.age(3,2,"rweibull", c(0.4,3))
plot(trialtheetrees[[1]])

##obtain 5 trees (with 40 tips) with taxa based and asymmetric 
#simulation, without extinction
mytree <-  sim.taxa(numbsim=5, n=40, m=40,  distributionspname="rweibull", distributionspparameters=c(0.1,1), distributionextname="rexp", distributionextparameters="0", symmetric = FALSE, complete=TRUE, labellivingsp="sp.", labelextinctsp="ext.", frac=1, sampling=2, gsa=FALSE)

#see sim.taxa and sim.age help pages for further details
}
