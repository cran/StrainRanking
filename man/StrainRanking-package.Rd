\name{StrainRanking-package}
\alias{StrainRanking-package}
\alias{StrainRanking}
\docType{package}
\title{
Package for ranking pathogen strains
}
\description{
Regression-based ranking of pathogen strains with respect to their contributions to natural epidemics, using demographic and genetic data sampled in the curse of the epidemics. 
}
\details{
\tabular{ll}{
Package: \tab StrainRanking\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-07-18\cr
License: \tab What license is it under?\cr
Depends: \tab methods\cr
}
To rank pathogen strains using the method of Soubeyrand et al. (2013), create a \link{DG object} (Demographic and Genetic data set) with one of the three construction functions (\link{DGobj.rawdata}, \link{DGobj.simul.regression} and \link{DGobj.simul.mechanistic}) and apply the \link{ranking.strains} function. Other construction functions returning a DG object might be written to extend the approach proposed by Soubeyrand et al. (2013). 
}
\author{
Soubeyrand, S., Tollenaere, C., Haon-Lasportes, E. and Laine, A.-L. 

Maintainer: \email{Samuel.Soubeyrand@avignon.inra.fr}
}
\references{
Soubeyrand, S., Tollenaere, C., Haon-Lasportes, E. and Laine, A.-L. (2013). Regression-based ranking of pathogen strains with respect to their contributions to natural epidemics. Research Report.
}

\keyword{ package }
\examples{
## Application of the ranking method to a real data set
data(powderymildew)
DGdata=DGobj.rawdata(demographic.coord=powderymildew$demographic.coord,
 genetic.coord=powderymildew$genetic.coord,
 demographic.measures=powderymildew$demographic.measures,
 genetic.frequencies=powderymildew$genetic.frequencies)
ranking.strains(DGobject=DGdata, bw=2, nb.mcsimul=10^3, plots=TRUE)

## Application of the ranking method to a data set simulated under the 
## mechanistic model
DGmech=DGobj.simul.mechanistic(sqrtn=10, size1=30, size2=10, theta=c(1.5,2,3), 
	beta=c(5,5), M=7, delta=0.2)
ranking.strains(DGobject=DGmech, bw=2, nb.mcsimul=10^3, plots=TRUE)
	
## Application of the ranking method to a data set simulated under the 
## regression model
DGreg=DGobj.simul.regression(sqrtn=10, size1=30, size2=10, theta=c(1.5,2,3), 
 alpha.function=generation.alpha.3strains, sigma=0.1)
ranking.strains(DGobject=DGreg, bw=2, nb.mcsimul=10^3, plots=TRUE)
}
