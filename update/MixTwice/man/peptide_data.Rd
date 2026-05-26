\name{peptide_data}
\alias{peptide_data}
\docType{data}
\title{Peptide array data example}
\description{
A high-density peptide microarray example to identify peptides for which antibody binding levels differ between control subjects and rheumatoid arthritis (RA) patients expressing a specific disease marker combination (i.e., CCP+RF+ RA). 
}
\usage{data("peptide_data")}
\format{
  A data frame with 152603 observations on the following 16 variables.
  \describe{
    The first 8 columns are RA patients and the remaining columns are from control subjects.
  }
}
\details{
Each row of the data (\code{rownames(peptide_data)}) is a probed length-12 peptide and each column of the data (\code{colnames(peptide_data)}) is a subject with distinct pseudo sample ID. The binding value is doubly-log transformed using natural base to stabilize variance.
}
\source{
Zheng, Zihao, et al. \emph{Disordered Antigens and Epitope Overlap Between Anti Citrullinated Protein Antibodies and Rheumatoid Factor in Rheumatoid Arthritis.} \bold{Arthritis & Rheumatology} 72.2 (2020): 262-272.
}
\references{
Zheng et al. \emph{MixTwice: Large scale hypothesis testing for peptide arrays by variance mixing}. Bioinformatics, 2021.
}
\examples{

#### load the RA data

data(peptide_data)

#### visualize the data

## each row is a peptide with unique peptide sequence

## each column is a subject with information on group and pseudo ID

colnames(peptide_data)

\donttest{
## z-score for peptide

get_zscore = function(x){
  
  n = length(x)
  
  t = t.test(x[1:(n/2)], x[(n/2 + 1):n], var.equal = TRUE)$statistic
  
  return(qnorm(pt(t, df = n-2)))
}

z = apply(peptide_data, 1, get_zscore)

## visualize the density of z-score

hist(z, probability = TRUE, 100, ylim = c(0,0.4), col = "blue")
lines(density(rnorm(10^5)), lwd =2)
}


}
\keyword{datasets}