\name{count.reads}
\alias{count.reads}
\title{Counting reads across intervals}
\description{Counts the number of reads starting at each position across given
genomic intervals}
\usage{count.reads(reads_fn, I, sb=NULL, binary=FALSE, sum.counts=FALSE)}
\arguments{
  \item{reads_fn}{filename of aligned reads in BAM format}
  \item{I}{a GRanges object giving valid genomic intervals}
  \item{sb}{a seqbias object}
  \item{binary}{if \code{TRUE}, return a 0-1 vector, otherwise return a vector
  counting the number of reads mapped to each position}
  \item{sum.counts}{if \code{TRUE} return the total read count for each interval}
}
\details{
  Given an indexed BAM file, this function counts the number of reads starting
  at each position of each provided interval. If a seqbias object is provided
  through the \code{sb} attribute, counts are corrected for sequence bias.
  The total read count for each interval is returned if \code{sum.counts} is
  \code{TRUE}.
  
  If \code{binary} is \code{TRUE} a 0-1 vector is returned instead, where
  positions at which no reads are mapped are 0, and those with one or more are
  1.
}
\value{
  A list of numeric vectors is returned, one for each interval provided. Each
  vector gives an integer count of the number of reads beginning on that
  position.
}
\note{
  The BAM file provided should be indexed with 'samtools index'.
}
\author{
    Daniel Jones
    \email{dcjones@cs.washington.edu}
}
\seealso{
    \code{\link{seqbias.predict}}
}
\examples{
  reads_fn <- system.file( "extra/example.bam", package = "seqbias" )

  I <- GRanges( c('seq1'), IRanges( c(1), c(5000) ), strand = c('-') )

  counts <- count.reads( reads_fn, I )
}


