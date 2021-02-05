p1  <- 1/5
p2  <- 1/4
p3  <- 1/3
p4  <- 13/60

AyIx  <- matrix(c(p1/(p1 + p3), p3/(p1 + p3), p2/(p2 + p4), p4/(p2 + p4)),
                nrow=2, byrow=TRUE)
AxIy  <- matrix(c(p1/(p1 + p2), p2/(p1 +p2), p3/(p3+p4), p4/(p3 + p4)),
                nrow=2, byrow=TRUE)

Xpi  <- 0                               #not needed
Ypi  <- 1                               #initial value for gibbsSeq

burnIn  <- 100
N  <- 1000
gibbsSeq  <- matrix(0, N, 2)

for (i in 1:N) {
    Xpi  <- rbinom(1,1, AxIy[2, Ypi + 1]) #recall indexing from 1
    Ypi  <- rbinom(1,1, AyIx[Xpi+1, 2])
    gibbsSeq[i, ]  <- c(Xpi, Ypi)
}

apply(gibbsSeq[(burnIn+1):N, ], 2, sum)/(N - burnIn)




                   
