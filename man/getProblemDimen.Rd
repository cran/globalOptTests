\name{getProblemDimen}
\alias{getProblemDimen}
\title{Get the length of the parameter vector expected by a given objective function.
}
\description{This function returns an integer representing the length of
  the parameter vector expected by a given objective function.  
}
\usage{
getProblemDimen(fnName)
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
\value{integer describing the length of the parameter vector expected by
  the objective function.  
}
\author{Katharine Mullen <katharine.mullen@stat.ucla.edu>}
\examples{
getProblemDimen("Ackleys")
getProblemDimen("AluffiPentini")
}
\keyword{optimization}