\name{broken.power.B}
\alias{broken.power.B}
\title{
ML estimation of broken-power law model with chosen number of breakpoints
}

\description{
This function conduct the maximum likelihood estimation of broken-power law model for astronomical sources. The methodology is described in \emph{Wong et. al., 2012}. This function requires the number of breakpoints (or pieces) to be chosen. When the number of pieces is set to 1, this function computes MLE for the ordinary \emph{log(N) - log(S)} model.
}

\usage{
broken.power.B(Y, A, B, taus0=NULL, betas0=NULL, b=rep(0,length(Y)), Nlim=30, Nsim=7000, Nburn=2000, thin=5, plot.check=FALSE, display=TRUE, NM.maxit=200, NM.funevals=300)
}

\arguments{
\item{Y}{a vector of observed counts}
\item{A}{a vector of values of effective areas}
\item{B}{number of pieces in the broken-power law model. If \code{B=1}, the function computes MLE for the ordinary \emph{log(N) - log(S)} model.}
\item{taus0}{a vector of initial values of breakpoints (including the starting point). If it is not provided, the function will generate the initial values automatically.}
\item{betas0}{a vector of initial values of `(negative) slopes'. If it is not provided, the function will generate the initial values automatically.}
\item{b}{a vector of expected blackground counts. If it is not provided, the function estimates the model with no background counts.}
\item{Nlim}{the maximum number of iterations for the interwoven EM algorithm}
\item{Nsim}{the number of simulations used in the Monte Carlo estimate of the Q function, stated in \emph{Wong, et. al.}}
\item{Nburn}{the number of burn-in used in the Monte Carlo estimate of the Q function, stated in \emph{Wong, et. al.}}
\item{thin}{the number of thining used in the Monte Carlo estimate of the Q function, stated in \emph{Wong, et. al.}}
\item{plot.check}{If \code{plot.check=TRUE}, some trace plots of flux are plotted for each iteration.}
\item{display}{If \code{display=TRUE}, progress of the grid search is displayed.}
\item{NM.maxit}{number of maximum iterations in the Nelder-Mead algorithm, used in the M step}
\item{NM.funevals}{number of maximum calls of (Q) function in the Nelder-Mead algorithm, used in the M step}
}

\value{
\item{pars}{a matrix of the estimated parameters for every iterations, after the sufficiency EM step; for internal use}
\item{pars1}{a matrix of the estimated parameters for every iterations, after the auxillary EM step; for internal use}
\item{sam.U}{for internal use}
\item{hat.pars}{for internal use}
\item{hat.taus}{the estimated values for taus}
\item{hat.betas}{teh estimated values for betas}
\item{EU}{for internal use}
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

