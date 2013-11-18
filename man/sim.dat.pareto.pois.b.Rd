\name{sim.dat.pareto.pois.b}
\alias{sim.dat.pareto.pois.b}
\title{
Simulate data from Broken-Power Law Model
}

\description{
This function simulates data from the broken-power law model.
}

\usage{
sim.dat.pareto.pois.b(n, taus, betas, A, b)
}

\arguments{
\item{n}{number of counts}
\item{taus}{a vector of breakpoints (including the starting point)}
\item{betas}{a vector of "negative slopes"}
\item{A}{a vector of values of effective areas}
\item{b}{a vector of expected blackground counts}
}

\value{
\item{Y}{the simulated counts}
\item{A}{the effective areas}
\item{S}{the simulated flux, usually not known by the observer}
\item{b}{the expected background counts}
}

\references{
Raymond K. W. Wong, Paul Baines, Alexander Aue, Thomas C. M. Lee and Vinay L. Kashyap (2012) \emph{}. Submitted.
}

\author{
Raymond K. W. Wong <raymondkww.dev@gmail.com>
}

\seealso{
\code{\link{broken.power}}
}

\examples{

}

