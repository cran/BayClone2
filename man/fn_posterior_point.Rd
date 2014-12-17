\name{fn_posterior_point}
\alias{fn_posterior_point}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fn_posterior_point function
%%  ~~function to do ... ~~
}
\description{This functions returns point estimates of the parameters such as a subclonal copy number matrix, a matrix of the number of copies with variant sequence in subclones, a matrix of composition weights of sampels in subclones and the expected read count with two copies (i.e posterior point estimate for L, Z, W, PHI, PI, P0 for a chosen value of C).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fn_posterior_point(CC, SS, TT, sam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{CC}{the number of subclones chosen by users. It should be less than or
equal to 10(including background subclone). This limitation is due to the permutaion.}
\item{SS}{the number of loci}
\item{TT}{the number of tissue samples}
\item{sam}{a list of MCMC samples returned from the function, BayClone2}
%%     ~~Describe \code{x} here~~
    
}
\details{The argument passed to this function, sam is a list returned from BayClone2; sam should be a list of posterior samples of random parameters (returned from the funtion, BayClone2); C, L, Z, w, th, phi, pi, p0_z, M and p
}

\value{This function returns 
	
C: the value of C passed to the function

L: a posteror point esitmate of L in a S*C matrix 

Z: a posteror point esitmate of Z in a S*C matrix 

w: a posteror point esitmate of w in a T*C matrix 

p0: a posteror point esitmate of p0 as a scalor

phi: a posteror point esitmate of phi in a vector of T

M: a posterior point estimate of M in a S*T matrix

p: a posterior point estimate of p in a S*T matrix
}
\references{
	J. Lee, P. Mueller, S. Sengupta, K. Gulukota, Y. Ji, Bayesian Inference for Tumor Subclones Accounting for Sequencing and Structural Variants (http://arxiv.org/abs/1409.7158)
}
\author{
	J. Lee (juheelee@soe.ucsc.edu) and S. Sengupta (subhajit06@gmail.com)
}

%\note{NA%TBD
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	export_N_n, BayClone2, fn_post_C
}
\examples{
##ILLUSTRATE BayClone2 WITH A SMALL SIMULATION.
###REPRODUCE SIMULATION 1 OF LEE ET AL.
library("BayClone2")

##READ IN DATA
data(BayClone2_Simulation1_mut)
data(BayClone2_Simulation1_tot)
##TOTAL NUMBER OF READS AT LOCUS s IN SAMPLE t
N <- as.matrix(BayClone2_Simulation1_tot)  
##NUMBER OF READS WITH VARIANT SEQUENCE AT LOCUS s IN SAMPLE t
n <- as.matrix(BayClone2_Simulation1_mut) 

S <- nrow(N)  # THE NUMBER OF LOCI (I.E. NUMBER OF ROWS OF N (AND n))
T <- ncol(N) #THE NUMBER OF TISSUE SAMPLES  (I.E. NUMBER OF COLUMNS OF N (AND n))

###################################
#HYPER-PARAMETER  ----SPECIFYING HYPERPARAMETER VALUES
######################################
#HYPER-PARAMETER
hyper <- NULL

#NUMBER OF SUBCLONES (GEOMETRIC DIST)
### C ~ GEOMETRIC(r) WHERE E(C)=1/r
hyper$r <- 0.2

#PRIOR FOR L
hyper$Q <- 3  #NUMBER OF COPIES -- q = 0, 1, 2, 3

##BETA-DIRICHLET
###PI_C | C ~ BETA-DIRICHLET (ALPHA/C, BETA, GAMMA)
hyper$alpha <- 2
hyper$beta <- 1
hyper$gam <- c(0.5, 0.5, 0.5)

#PRIOR FOR PHI--TOTAL NUMBER OF READS IN SAMPLE T
###PHI_T ~ GAMMA(A, B)
hyper$b <- 3
hyper$a <- median(N)*hyper$b

#PRIOR FOR P_O
###P0 ~ BETA(a, b)
hyper$a_z0 <- 0.3
hyper$b_z0 <- 5

#PRIOR FOR W
##W_T | L ~ DIRICHLET(D0, D, ..., D) WHERE W_T=(w_t0, w_t1, ..., w_tC)
hyper$d0 <- 0.5
hyper$d <- 1

#WE USE THE MCMC SIMULATION STRATEGY PROPOSED IN LEE AT EL (2014)
n.sam <- 10000;  ##NUMBER OF SAMPLES THAT WILL BE USED FOR INFERENCE
##NUMBER OF SAMPLES FOR BURN-IN 
#(USE THIS FOR A TRAINING DATA---FOR DETIALS, SEE THE REFERENCE)
burn.in <- 6000  

##############################################
###WE CONSIDER C BETWEEN 1 AND 15 IN ADDITION TO BACKGROUND SUBCLONE
####Max_C AND Min_C SPECIFIES VALUES OF C FOR POSTERIOR EXPLORATION
Min_C <- 2  ##INCLUDING THE BACKGROUND SUBCLONE
Max_C <- 16  ##INCLUDING THE BACKGROUND SUBCLONE


#################################################################
##DO MCMC SAMPLING FROM BAYCLONE2!
#################################################################
##THE LAST ARGUMENT (0.025) IS THE MEAN PROPORTION FOR THE TRAINING DATASET (SPECIFIED BY USERS)
##IT WILL BE USED TO SPLIT INTO TRAINING AND TEST DATASETS
##FOR DETAILS, SEE THE REFERENCE LEE AT EL (2014)
##TO RUN, COMMENT IN THE LINE BELOW (WARNING! THIS MAY TAKE APPROXIMATELY 30 MINUTES)
#set.seed(11615)
#MCMC.sam <- BayClone2(Min_C, Max_C, S, T, burn.in, n.sam, N, n, hyper, 0.025)


#################################################################
#COMPUTE THE POSTERIOR MARGINAL DIST OF C (THE NUMBER OF SUBCLONES)
#################################################################
##TO RUN, COMMENT IN THE LINE BELOW
#post_dist_C <- fn_post_C(MCMC.sam$C, Min_C, Max_C)

######################################################################################
####WE FIND POSTERIOR POINT ESTIMATES OF L, Z, W, PHI, PI, P0 FOR A CHOSEN VALUE OF C
######################################################################################
##THE FIRST ARGUMENT (3) IS A VALUE OF C CHOSEN BY USERS
#C IS THE NUMBER OF SUBCLONES INCLUDING THE BACKGROUPD SUBCLONE
##THE CHOSE VALUE OF C SHOULD BE LESS THAN OR EQUAL TO 10 (INCLUDING THE BACKGROUND SUBCLONE)
#DUE TO THE PERMUTATION (FOR DETAILS, SEE SEE THE REFERENCE LEE AT EL (2014))
##TO RUN, COMMENT IN THE LINE BELOW (WARNING! THIS MAY TAKE ARPPOXIMATELY 15 MINUTES)
#point.est <- fn_posterior_point(3, S, T, MCMC.sam)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Summary-of-Posterior-Samples}% __ONLY ONE__ keyword per line
