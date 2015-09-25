\name{str.seir}

\alias{str.seir}

\title{Describe the Structure of a seir Object.}

\description{Describe the structure of a seir object.}

\usage{str(object, ...)}

\arguments{
  \item{object}{a seir object.}
  \item{\dots}{further arguments.}
}

\details{Compactly display the internal structure of an object, of the 'seir' class. It displays the names and the values of all parameters for each subpopulation in a metapopulation.}

\value{'str.seir' does not return anything, for efficiency reasons.  The obvious side effect is output to the terminal.}

\author{TRAN Thi Cam Giang}

\seealso{'str.seir' function in the 'dizzys' package.}

\examples{
	obj<-seir()
	str(obj)
}

\keyword{ seir model }
\keyword{ R object }
