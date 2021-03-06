% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BCC1997.R
\name{BCC}
\alias{BCC}
\title{Calculation of Option Prices Based on a Universal Solution}
\usage{
BCC(kappav, kappar, thetav, thetar, sigmav, sigmar, muj, sigmaj, rho, lambda,
  S0, K, V0, R0, t)
}
\arguments{
\item{kappav}{Speed of convergence on variance}

\item{kappar}{Speed of convergence on risk free rate}

\item{thetav}{Long-term variance}

\item{thetar}{Long-term risk free rate}

\item{sigmav}{Volatility of variance}

\item{sigmar}{Volatility of risk free rate}

\item{muj}{Jump size}

\item{sigmaj}{Volatility of jumps}

\item{rho}{Correlation between underlying price and variance}

\item{lambda}{Jump intensity}

\item{S0}{Initial/Current underlying price}

\item{K}{Strike price}

\item{V0}{Initial/Current variance}

\item{R0}{Initial/Current risk free rate}

\item{t}{Time to maturity}
}
\value{
Call: return the price of the European call oprion

Put: return the price of the European put oprion
}
\description{
This is a function to calculate the prices of European options based on the universal solution provided by Bakshi, Cao and Chen (1997) <doi:10.1111/j.1540-6261.1997.tb02749.x>. This solution takes stochastic volatility, stochastic interest and random jumps into consideration. Please cite their work if this package is used.
}
\note{
Please notice each parameter has its "reasonable range". e.g. volatilities cannot be zero or smaller than zero, please input 0.0000001 when they are zero.
}
\examples{
BCC(kappav=0,kappar=0,thetav=0,thetar=0,sigmav=0.0000001,sigmar=0.0000001,muj=0,
         sigmaj=0.0000001,rho=0,lambda=0,S0=100,K=100,V0=0.04,R0=0.01,t=1)
BCC(kappav=0.5,kappar=0,thetav=0.025,thetar=0,sigmav=0.09,sigmar=0.0000001,muj=0,
         sigmaj=0.0000001,rho=0.1,lambda=0,S0=100,K=100,V0=0.04,R0=0.01,t=1)
}
