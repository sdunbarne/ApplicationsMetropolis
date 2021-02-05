makePicture  <- function(N) {
    ## Make an original N x N  picture of 0s and 1s to unblur,
    ## works nicely if N is multiple of 40, for tests, multiple of 8   
    picture  <- matrix(0, N, N)

    spot1Top  <- floor( N/2 )
    spot1Bottom  <- floor( (7/8) * N )
    spot1Left  <- floor( (1/8) * N )
    spot1Right  <- floor( (1/4) * N )

    spot2Top  <- floor( (1/8) * N )
    spot2Bottom  <- floor( (1/4) * N )
    spot2Left  <- floor( (3/4) * N )
    spot2Right  <- floor( (7/8) * N )

    picture[spot1Top:spot1Bottom, spot1Left:spot1Right]  <-  1
    picture[spot2Top:spot2Bottom, spot2Left:spot2Right]  <-  1

    picture
}

makeNoise  <- function(N, sigma=1) {
    ## N x N matrix of Gaussian noise
    noise  <- matrix( rnorm( N*N, mean=0, sd=sigma), N, N)
}

roundPicture  <- function(M) {
    ## Since the values of \( \omega^{\text{blurred}} \) are real
    ## numbers, the resulting image is determined by rounding each
    ## value to  the nearest value in \( S = {0,1 \).
    ifelse( M < 1/2, 0, 1)
}

makeBlurredPicture  <- function(N, sigma=0.3) {
    ## make a speckled, blurred image to work with
    ## default sigma  = 0.3 is from experimentation, gives good pic
    blurPict <- roundPicture( makePicture(N) + makeNoise(N, sigma=0.3))
}

surroundNeighbors  <- function(site, N) {
    ## Note that these are periodic neightbors!, may want to change
    belowSite  <- if (site[1] == N) c(1, site[2]) else c(site[1]+1, site[2])
    leftSite  <- if (site[2] == 1) c(site[1], N) else c(site[1], site[2]-1)
    aboveSite  <- if (site[1] == 1) c(N, site[2]) else c(site[1]-1, site[2])
    rightSite  <- if (site[2] == N) c(site[1], 1) else c(site[1], site[2]+1)
    neighbors  <- c(belowSite, leftSite, aboveSite, rightSite)
}

probSite  <- function(site, pict, blurred, N, sigma=0.3) {
    ij  <-  surroundNeighbors(site, N)
    
    p1  <- exp( ( -1/(2 * sigma^2) *
                 (1 - blurred[site[1], site[2]] )^2
                ) -
                 1 * ( pict[ ij[1], ij[2] ] +
                       pict[ ij[3], ij[4] ] +
                       pict[ ij[5], ij[6] ] +
                       pict[ ij[7], ij[8] ]
                      )
               )
    p0  <- exp( ( -1/(2 * sigma^2) *
                 (0 - blurred[site[1], site[2]] )^2
                ) -
                 0 * ( pict[ ij[1], ij[2] ] +
                       pict[ ij[3], ij[4] ] +
                       pict[ ij[5], ij[6] ] +
                       pict[ ij[7], ij[8] ]
                      )
               )


    pr  <-  p0/(p0 + p1)
}

gibbsSample  <-  function(site, pict, blurred, sigma=0.3) {
    ## The Gibbs sampling procedure, pick a pixel value to
    ## balance straying too far from data (blurred), but
    ## also trying to stay aligned with neighbors, as in Ising
    x  <- runif(1)
    p   <- probSite(site, pict, blurred, N, sigma=0.3)
    if ( x < p ) {
        k  <- 0
    } else {
        k  <- 1
    }
    k
}

rasterScan  <- function(pict, blurred, sigma=0.3) {
    NR = NROW(pict)
    NC = NCOL(pict)
    
    for (i in 1:NR) {
        for (j in 1:NC) {
            site  <-  c(i,j)
            pict[i, j]  <- gibbsSample(site, pict, blurred, sigma=0.3)
        }
    }
    pict
}

N  <- 16
nIter  <-  5000

testPict  <- makePicture(N)
testblurredPict  <-  makeBlurredPicture(N)

pict  <- testblurredPict

for (m in 1:nIter) {
    newpict   <- rasterScan(pict, testblurredPict)
    pict  <- newpict
}

par(mfrow=c(1,2))
image(testPict)
image(pict)

## NAME: gibbsImage.R
## USAGE: within R, at interactive prompt
##        source("gibbsImage.R")
## REQUIRED ARGUMENTS: none
## OPTIONS:  none
## DESCRIPTION: Uses Gibbs sampling to attempt to restore a
##              speckled and blurred image
## DIAGNOSTICS: none
## CONFIGURATION AND ENVIRONMENT: base R
## DEPENDENCIES: base R
## INCOMPATIBILITIES: none known
## PROVENANCE: Steve Dunbar
## BUGS AND LIMITATIONS: Slow, doesn't seem to do a good restoration.
## FEATURES AND POTENTIAL IMPROVEMENTS: slow, needs speed up.
## AUTHOR:  Steve Dunbar
## VERSION: Version 1.0 as of Tue 05 Jan 2021 08:44:17 AM CST
## KEYWORDS: Gibbs sampling, image restoration, Metropolis alooirhtm

