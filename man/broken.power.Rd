\name{broken.power}
\alias{broken.power}
\title{
ML estimation of broken-power law model with automatic selection of number of breakpoints
}

\description{
This function conduct the maximum likelihood estimation of broken-power law model for astronomical sources. The methodology is described in \emph{Wong et. al., 2012}. This function does not require the number of breakpoints (or pieces) to be chosen.
}

\usage{
broken.power(Y, A, maxB=3, b=rep(0,length(Y)), Nlim=30, Nsim1=7000, Nburn1=2000, Nsim2=10000, Nburn2=2000, thin1=5, thin2=5, plot.check=FALSE, display=FALSE, PPmulti=TRUE, cpus=4, NM.maxit=200, NM.funevals=300, details=FALSE,useC=TRUE)
}

\arguments{
\item{Y}{a vector of observed counts}
\item{A}{a vector of values of effective areas}
\item{maxB}{maximum number of pieces in the broken-power law model the function will compare. If \emph{maxB=3}, the function compares model with B=1, 2 and 3 and then choose the best fit model.}
\item{b}{a vector of expected blackground counts. If it is not provided, the function estimates the model with no background counts.}
\item{Nlim}{the maximum number of iterations for the interwoven EM algorithm}
\item{Nsim1}{the number of simulations used in the Monte Carlo estimate of the Q function for the interwoven EM algorithm, stated in \emph{Wong, et. al.}}
\item{Nsim2}{the number of simulations used in the Monte Carlo estimate of the log likelihood for the power posterior method, stated in \emph{Wong, et. al.}}
\item{Nburn1}{the number of burn-in used in the Monte Carlo estimate of the Q function for the interwoven EM algorithm, stated in \emph{Wong, et. al.}}
\item{Nburn2}{the number of burn-in used in the Monte Carlo estimate of the log likelihood for the power posterior method, stated in \emph{Wong, et. al.}}
\item{thin1}{the number of thining used in the Monte Carlo estimate of the Q function for the interwoven EM algorithm, stated in \emph{Wong, et. al.}}
\item{thin2}{the number of thining used in the Monte Carlo estimate of the log likelihood for the power posterior method, stated in \emph{Wong, et. al.}}
\item{plot.check}{If \code{plot.check=T}, some trace plots of flux are plotted for each iteration.}
\item{display}{If \code{display=TRUE}, progress of the grid search is displayed.}
\item{PPmulti}{If \code{PPmulti=TRUE}, parallel computing will be applied for the computation of the log likelihood through power posterior method.}
\item{cpus}{number of cores used for the parallel computing; only useful when \code{PPmulti=T}}
\item{NM.maxit}{number of maximum iterations in the Nelder-Mead algorithm, used in the M step}
\item{NM.funevals}{number of maximum calls of (Q) function in the Nelder-Mead algorithm, used in the M step}
\item{details}{If \code{details=TRUE}, the function returns a more detailed output.}
\item{useC}{If \code{useC=TRUE}, C code is used, otherwise R code is used.}
}

\value{
\item{B}{estimated number of pieces (by BIC)}
\item{hat.taus}{estimated taus (by BIC)}
\item{hat.betas}{estimated betas (by BIC)}
\item{est.pars}{a list containing all estimates of taus and betas for different B}
\item{IC}{a matric containing the AIC (column 1) and BIC (column 2) for different B}
\item{bestB}{a vector of numbers of pieces selected by AIC (first entry) and BIC (second entry)}
\item{loglike}{a vector of estimated likelihood for different B}
\item{out}{If \code{details=TRUE}, this contains a list of the output from \code{\link{broken.power.B}} for different B}
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
# TODO
}

