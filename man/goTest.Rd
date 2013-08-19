\name{goTest}
\alias{goTest}
\title{Call an objective function
}
\description{This function is used to access one of 50 objective
  functions that may be useful for benchmarking performance of global
  optimization algorithms. 
}
\usage{
goTest(par, fnName, checkDim = TRUE)
}
\arguments{
  \item{par}{vector of reals representing the parameter vector at which
    to evaluate the objective function
}
  \item{fnName}{a character vector representing the name of the
    objective function to use.  Options are elements of 
 \preformatted{c("Ackleys", "AluffiPentini", "BeckerLago",
      "Bohachevsky1", "Bohachevsky2", "Branin",
      "Camel3", "Camel6", "CosMix2", "CosMix4",
      "DekkersAarts", "Easom", "EMichalewicz",
      "Expo", "GoldPrice", "Griewank", "Gulf",
      "Hartman3", "Hartman6", "Hosaki", "Kowalik",
      "LM1", "LM2n10", "LM2n5", "McCormic",
      "MeyerRoth", "MieleCantrell", "Modlangerman",
      "ModRosenbrock", "MultiGauss", "Neumaier2",
      "Neumaier3", "Paviani", "Periodic",
      "PowellQ", "PriceTransistor", "Rastrigin",
      "Rosenbrock", "Salomon", "Schaffer1",
      "Schaffer2", "Schubert", "Schwefel",
      "Shekel10", "Shekel5", "Shekel7",
      "Shekelfox5", "Wood", "Zeldasine10",
      "Zeldasine20")}.
      Unique partial matches to these elements
      are accepted. 
}   
\item{checkDim}{a boolean value that is \code{TRUE} if the length of the
  parameter vector should be checked to make sure it matches the
  expected value.  If set to \code{FALSE}, the function runs slightly faster.
}}
\value{a real that representes the value of the objective function
  \code{fnName} at the parameter vector \code{par}. 

}
\references{
  \code{C} code for the functions implemented here was originally by 
\cite{
M. Montaz Ali, Charoenchai Khompatraporn and Zelda B. Zabinsky, Journal
of Global Optimization, 2005, 31, 4, 635-672. 
} 
The \code{C} implementations were downloaded in a form that was modified from
the paper above, and modified further for incorporation into this
package from a site maintained by GAMS Development Corp. and GAMS
Software GmbH:
\url{http://www.gamsworld.org/performance/selconglobal/selcongloballib.htm}.

Note that the factor of 0.02 in Ackley's function has been changed to
  the value commonly found in the literature, 0.2.  Also, Storn's
Tchebychev Problem in 9 and 17 dimensions was not included, since the
global minima of the imlementation of these functions does not
correspond to the global minima reported in the above paper.
}

\author{Katharine Mullen \email{katharine.mullen@stat.ucla.edu}}

\examples{
goTest(fnName="Ackleys", par=rnorm(10))
goTest(fnName="AluffiPentini", par=c(1,2))
goTest(fnName="AluffiPentini",
par=rep(1,getProblemDimen("AluffiPentini")))

## use in an optimization via 'optim'
optim(par=c(1,2), fn=goTest, fnName="AluffiPentini")
}
\keyword{optimization}