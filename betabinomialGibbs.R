alphaParam  <- 2
betaParam  <- 4
n  <- 16
trials  <- 500

Xpi  <- rep(8, trials)                  #not needed
Ypi  <- rep(0.5, trials)                #initial value for gibbsSeq

gibbsIter  <- 10

for (i in 1:gibbsIter) {
    Xpi  <- rbinom(trials, 16, Ypi)     #note vectorization!
    Ypi  <- rbeta(trials, Xpi + alphaParam, n - Xpi + betaParam)
}

hist(Xpi, breaks=n)

## NAME: betabinomialGibbs.R
##       Use Gibbs sampler to generate samples from beta-binomial pdf
## USAGE: within R, at interactive prompt
##        source("betabinomial.R")
## REQUIRED ARGUMENTS: none
## OPTIONS: none
## DESCRIPTION: 
## Use Gibbs sampler to generate samples from beta-binomial pdf
## and then plot a histogram of the samples             
## DIAGNOSTICS: none
## CONFIGURATION AND ENVIRONMENT: base R
## DEPENDENCIES: base R
## INCOMPATIBILITIES: none known
## PROVENANCE: Steve Dunbar, December 31, 2020
## BUGS AND LIMITATIONS: none known
## FEATURES AND POTENTIAL IMPROVEMENTS: excellent candidate for 
##        (naive) parallelization
## AUTHOR:  Steve Dunbar
## VERSION: Version 1.0 as of December 1, 2020
## KEYWORDS: Gibbs sampler, binomial, beta, beta-binomial, gamma





                   
