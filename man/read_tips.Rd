\name{read_tips}
\alias{read_tips}
\title{Read the tip states and labels from a file}
\usage{
read_tips(file, character=1, sep="\t")
}
\arguments{
    \item{file}{The path to the file containing the data.}
    \item{character}{If the file contains more than one character, you can specify which one to read by the column index.}
    \item{sep}{The csv separator, defaults to TAB.}

}
\description{
    Read a CSV file given as argument and return a matrix representing the states of the tips and its labels (taxa names).

    The CSV file must be in a specific format. The first column contain the taxa names, each one in a different row. Every subsequent column is a character, and the state of the character is specified by each taxa on each row, forming a matrix.

    You can have as many characters (columns) as you want, but you must specify one to be read, using the \code{character} argument, which is an integer representing the column of the character (first character is in column 1).

    When some taxon have an ambiguous state, put all states side-by-side, without an separator. For instance, if taxa X can be at both states 'a' and 'b', write 'ab' in the file.
}
\details{
    If you decide to use spaces to separate taxa names from characters make sure the taxa names doesn't have spaces. For instance: "A_grahamae_Agra_Cur" is correct, "A grahamae Agra Cur" is wrong.
}
\value{
    {A matrix with the taxa label as the row names and a number of columns equal to the number of possible states of the character. The values of the matrix can be one, when the taxon can be at the state, and 0 otherwize.}
}

%\references{
%    References
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
%\seealso{
%    See Also
%}
%\keyword{phylogenetics}
\keyword{utilities}
