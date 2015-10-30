\name{summary.seir}

\alias{summary.seir}


\title{Object Summaries}

\description{generic function used to display the internal structure of an object of 'seir' class.}

\usage{summary.seir(object, ...)}

\arguments{
  \item{object}{an objecct of the 'seir' class.}
  \item{\dots}{further arguments.}
}

\details{Display in detail the internal structure of an object of the 'seir' class. It displays the names, the values, the signification of all parameters for each subpopulation in a metapopulation.
}

\value{'str.seir' does not return anything, for efficiency reasons.  The obvious side effect is output to the terminal.}

\author{Tran Thi Cam Giang.}

\seealso{'summary.seir' function of the 'dizzys' package.}

\examples{
	obj<-seir()
	summary(obj)
}

\keyword{ summary }
\keyword{ seir class }
