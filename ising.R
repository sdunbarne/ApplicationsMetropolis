downConfig <- function(N) { matrix(-1,N,N) } # 1 of 2 lowest energy configs
upConfig <- function(N) { matrix(1,N,N) }    # other lowest energy configs

randomConfig <- function(N) {
  config <- matrix(2*sample(0:1, N*N, replace=TRUE)-1, N,N)
}

randomSite  <- function(N) {
    rns  <- runif(2)
    i <- floor( N* rns[1] ) + 1
    j <- floor( N* rns[2] ) + 1
    c(i,j)
}
## This is a fairly efficient way to make a random site, but
## this is the hotspot in the script, occupying most of the runtime.

boltzWeight <- function(nrg, kTemp) { exp(- nrg/kTemp) }
## This has the negative sign, same as both Richey and Schlusser
## But this does not use J, already in energy, consistent with Richey (but not Schlusser)
## Uses kTemp, Boltzmann constant times temperature, with units of energy
## so the energy units cancel.

surroundNeighbors  <- function(site) {
    belowSite  <- if (site[1] == N) c(1, site[2]) else c(site[1]+1, site[2])
    leftSite  <- if (site[2] == 1) c(site[1], N) else c(site[1], site[2]-1)
    aboveSite  <- if (site[1] == 1) c(N, site[2]) else c(site[1]-1, site[2])
    rightSite  <- if (site[2] == N) c(site[1], 1) else c(site[1], site[2]+1)
    neighbors  <- c(belowSite, leftSite, aboveSite, rightSite)
}

metropolisMCStep  <- function(config, N) {
    site  <- randomSite(N)
    configAtSite  <- config[ site[1], site[2] ]
    configPrimeAtSite  <- configAtSite * (-1)

    neigh  <- surroundNeighbors(site) 
    deltaEnergy <- -2*J * (configPrimeAtSite - configAtSite)*
        (config[ neigh[1], neigh[2] ] +
         config[ neigh[3], neigh[4] ] +
         config[ neigh[5], neigh[6] ]+
         config[ neigh[7], neigh[8] ]
        )

    ## IF deltaEnergy < 0, so energy(configPrime) < energy(config)
    ## so p = boltzWeight(deltaEnergy, kTemperature) > 1,
    ## THEN  runif(1) is certain to be less than p so certainly
    ## change config to configPrime
    ## ELSE deltaEnergy > 0, so 
    ## energy(configPrime) > energy(config), so 
    ## p = boltzWeight(deltaEnergy, kTemperature) < 1,
    ## so change only on chance that runif(1) < p
    p  <- boltzWeight(deltaEnergy, kTemperature)
    
    if ( runif(1) < p ) {
        config[ site[1], site[2] ]  <- configPrimeAtSite
        config
    } else {
        config
    }
}

J  <- 1

kTemperature  <- 1.7

N  <- 40
mcmcSteps  <- 40000

config  <- randomConfig(N)

for (i in 1:mcmcSteps ) {
    config  <- metropolisMCStep(config, N)
}

image(1:N, 1:N, config[N:1, ], 
      col=gray.colors(2),
      xlab="", ylab="", xaxt="n", yaxt="n")
## config[N:1, ] so that the plot corresponds to matrix order
## i.e. first row of matrix is at top of plot, last row at bottom
