\name{TreeSimGM-package}
\alias{TreeSimGM-package}
\alias{TreeSimGM}
\docType{package}
\title{Simulating Phylogenetic Trees under General Bellman Harris and Lineage Shift Model}
\description{
Provides a flexible simulation tool for phylogenetic trees under a general model for speciation and extinction. Trees with a user-specified number of extant tips, or a user-specified stem age are simulated. It is possible to assume any probability distribution for the waiting time until speciation and extinction. Furthermore, the waiting times to speciation / extinction may be scaled in different parts of the tree, meaning we can simulate trees with clade-dependent diversification processes. At a speciation event, one species splits into two. We allow for two different modes at these splits: (i) symmetric, where for every speciation event new waiting times until speciation and extinction are drawn for both daughter lineages; and (ii) asymmetric, where a speciation event results in one species with new waiting times, and another that carries the extinction time and age of its ancestor. The symmetric mode can be seen as an vicariant or allopatric process where divided populations suffer equal evolutionary forces while the asymmetric mode could be seen as a peripatric speciation where a mother lineage continues to exist.
}
\details{
\tabular{ll}{
Package: \tab TreeSimGM\cr
Type: \tab Package\cr
Version: \tab 2.5\cr
Date: \tab 2020-03-04\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
}
TreeSimGM returns a list of phylogenetic trees. The trees can be simulated in different ways: (i) by age 'sim.age' (i.e. the simulation is started with one species and the user specifies a time, which is the stem age after which the simulation should stop); (ii) by taxa 'sim.taxa' (i.e. the user specifies the number of extant tips in the final tree). The user always specifies a distribution and parameters for the waiting time until speciation / extinction. Moreover, the user needs to specify the mode (i.e. symmetric or asymmetric), thus, there are four major possible combinations: (a) sim.age with parameter symmetric = TRUE; (b) sim.age with parameter symmetric = FALSE (i.e. asymmetric); (c) sim.taxa with parameter symmetric = TRUE; (d) sim.taxa with parameter symmetric = FALSE (i.e. asymmetric). Additionally, the user can choose if lineage-specific changes shall occur, through scaling in the waiting times to speciation or extinction events. The scaling can happen for (1) speciation, (2) extinction or (3) both. Setting the probability for a new species to have a scaled speciation or extinction times ('shiftspprob' and 'shiftextprob') to a value between ]0,1] will automatically activate the lineage-specific change model. If a change happens, a value is taken from the speciation / extinction waiting time distribution and is multiplied by a stochastic scaling factor. By default, both speciation and extinction shift frequencies are set to zero, meaning that if not stated otherwise, no shifts are considered.
}
\author{
Oskar Hagen, Tanja Stadler

Maintainer: Oskar Hagen <oskar@hagen.bio>
}
\references{

O. Hagen and T. Stadler (2017). TreeSimGM: Simulating phylogenetic trees under general Bellman Harris models with lineage-specific shifts of speciation and extinction in R. Methods in Ecology and Evolution.

O. Hagen, K. Hartmann, M. Steel and T. Stadler (2015). Age-dependent Speciation can Explain the Shape of Empirical Trees. Systematic Biology. 64, v.3, p.432-440.
}

\concept{tree}
\concept{phylogenies}
\concept{phylogenetic trees}
\concept{simulation}
\concept{general model} 
\concept{macroevolution}

\seealso{
\code{\link[TreeSim]{TreeSim}},
\code{\link[ape]{ape}}, \code{sim.age}, \code{sim.taxa}, \code{track.shift}
}
\examples{
##plots the first tree of a list of two simulated trees
#simulation based on age, this case = 3. 
#Note that by default, symmetric = TRUE
#i.e. we simulate under symmetric speciation and under no extinction.
library("TreeSimGM")
trialtrees <- sim.age(3,4,"rweibull(0.4,3)")
plot(trialtrees[[1]])

##obtain 5 trees with 40 tips under asymmetric speciation and with extinction
#simulation, without extinction
mytree <-  sim.taxa(numbsim=5, n=40, m=40,  waitsp="rweibull(0.1,1)", 
waitext="rexp(0.5)", symmetric = FALSE)

#see sim.age and sim.taxa help pages for further details
}
