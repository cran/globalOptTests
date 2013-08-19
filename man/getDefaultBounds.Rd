\name{getDefaultBounds}
\alias{getDefaultBounds}
\title{Get the default box constraints for the objective function.}
\description{This function returns real vectors representing the default
  lower and upper bounds (box constraints) for the given objective
  function. 
}
\usage{
getDefaultBounds(fnName)
}
\arguments{
  \item{fnName}{a character vector naming the objective function to
    call. Options are elements of
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
}
\value{a list containing elements
\item{lower}{vector of reals representing lower bounds}
\item{upper}{vector of reals representing upper bounds}
}
\author{Katharine Mullen <katharine.mullen@stat.ucla.edu>}
\examples{
getProblemDimen("Ackleys")
getProblemDimen("AluffiPentini")
}
\keyword{optimization}
