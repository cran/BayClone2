\name{export_N_n}
\alias{export_N_n}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{export_N_n function
%%  ~~function to do ... ~~
}
\description{This function takes a VCF file as input and generates
two output files (names are specified by users) which contains the total 
number of reads and the number of reads that bear a mutated sequence, respectively, 
at a particular locus in a specific tissue sample.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
export_N_n(vcf_in_file,N_tot_out_file,n_alt_out_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{vcf_in_file}{the name of the input VCF file.}
\item{N_tot_out_file}{the output file name for total reads count (character string)}
\item{n_alt_out_file}{the output file name for mutant reads count (character string)}
%%     ~~Describe \code{x} here~~
    
}
\details{
Input VCF filename is passed as the first argument.
Second and third arguments are the name of the output files.

Second file contains total number of reads at loci in all samples.

Third file contains total number of mutant reads at loci in all samples. 

NOTE: Each row in both the output files represents one particular loci for all the samples.  The output files may be used as data, N and n, for the function, BayClone2.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The function generates two output files.
	
  -N_tot_out_file: the output file that contains a (S x T) matrix, where S is the number of loci and T is the number of samples. 
Each element represents total number of reads at a locus in a specific tissue sample.

  -n_alt_out_file: the output file that contains
a (S x T) matrix, where S is the number of loci and T is the number of samples. Each element represents number of mutated reads at a locus in a tissue sample.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
	J. Lee, P. Mueller, S. Sengupta, K. Gulukota, Y. Ji, Bayesian Inference for Tumor Subclones Accounting for Sequencing and Structural Variants (http://arxiv.org/abs/1409.7158)	
	
	Sengupta S, Gulukota K, Lee J, Mueller, P, Y. Ji, BayClone: Bayesian Nonparametric Inference of Tumor Subclones Using NGS Data. Conference paper accepted for PSB 2015 and oral presentation
}
\author{
	J. Lee (juheelee@soe.ucsc.edu) and S. Sengupta (subhajit06@gmail.com)
}
%\note{NA%TBD
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
		BayClone2, fn_post_C, fn_posterior_point
}
\examples{
### Illustrate the functionality of the function export_N_n with an example
# please put the appropriate file names to run this function.
#library("BayClone2")
##INPUT FILE: test_Data.VCF
##OUTPUT FILES: N_tot.txt and n_alt.txt
#export_N_n("test_Data.VCF","N_tot.txt","n_alt.txt")

##LOAD THE OUTPUT FILES -- THE OUTPUT FILES CAN BE USED AS DATA FOR BAYCLONE2.
#N <- read.table("N_tot.txt")
#n <- read.table("n_alt.txt")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{VCF}
