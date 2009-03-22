# This file contains to sets of functions to calculate the mean and variance of D2
#
# The first set of functions (by Conrad Burden) is for the general case with
# mismatches, but only work for symetric Berouilli DNA, and does not include
# the off-diagonal part ofthe accordion
#
# The second set computes the exact variance (including off-diagonal accordion)
# for any iid text but only works for exact word matches


# The following functions compute the variance and mean of D2 for approximate
# word matches. Code by Conrad Burden.


g <- function(t,k,eta,c)
{
    h <- (1/4^k)*(1 - eta)^c*(1 + eta)^(k - c)
    iii <- 0:(k-t)
    v <- ((3 + eta)/(1 - eta))^(c - iii) * ((3 - eta)/(1 + eta))^(t + iii - c)
    u <- sum(choose(c,iii)*choose(k - c,k - t - iii)*v)
    g <- h*u
    g
}

G <- function(t,k,eta,cc)
{
    gg <- c()
    for (j in 0:t)
    {
        gg <- c(gg, g(j, k, eta, cc))
    }
    G <- sum(gg)
    G
}

mean.D2 <- function(t,k,eta,nA,nB)
{
    mean.D2 <- nA*nB*pbinom(t,k,(3 - eta^2)/4)
    mean.D2
}

fr.of.W <- function(r,cc,eta,t,k)
{
    fr.of.W <- 0
    for (l in 0:min(r,t))
    {
        fr.of.W <- fr.of.W + dbinom(l,r,(3 - eta^2)/4)*G(t - l,k - r,eta,cc)
    }
    fr.of.W
}

E.of.fr.of.W <- function(r,eta,t,k)
{
    E.of.fr.of.W <- 0
    for (cc in 0:(k - r))
    {
        E.of.fr.of.W <- E.of.fr.of.W + dbinom(cc,k - r,(1 - eta)/2)*fr.of.W(r,cc,eta,t,k)
    }
    E.of.fr.of.W
}

E.of.fr.of.W.sq <- function(r,eta,t,k)
{
    E.of.fr.of.W.sq <- 0
    for (cc in 0:(k - r))
    {
        E.of.fr.of.W.sq <- E.of.fr.of.W.sq + dbinom(cc,k - r,(1 - eta)/2)*fr.of.W(r,cc,eta,t,k)^2
    }
    E.of.fr.of.W.sq
}

Var.of.fr.of.W <- function(r,eta,t,k)
{
    Var.of.fr.of.W <- E.of.fr.of.W.sq(r,eta,t,k) - E.of.fr.of.W(r,eta,t,k)^2
    Var.of.fr.of.W
}

var.D2.crabgrass <- function(t,k,eta,nA,nB)
{
    var.D2.crabgrass <- Var.of.fr.of.W(0,eta,t,k)
    for (r in 1:(k - 1))
    {
        var.D2.crabgrass <- var.D2.crabgrass + 2*Var.of.fr.of.W(r,eta,t,k)
    }
    var.D2.crabgrass <- nA*nB*(nA + nB - 4*k + 2)*var.D2.crabgrass    
    var.D2.crabgrass
}

cov.Yu.Yv.acc.diag <- function(r,eta,t,k)
{
    top <- min(k - r,t)
    sum(dbinom(0:top,k - r,(3 - eta^2)/4) * pbinom(t - (0:top),r,(3 - eta^2)/4)^2) - pbinom(t,k,(3 - eta^2)/4)^2
}

var.D2.acc.diag <- function(t,k,eta,nA,nB)
{
    var.D2.acc.diag <- cov.Yu.Yv.acc.diag(0,eta,t,k)
    for (r in 1:(k - 1))
    {
        var.D2.acc.diag <- var.D2.acc.diag + 2*cov.Yu.Yv.acc.diag(r,eta,t,k)
    }
    var.D2.acc.diag <- nA*nB*var.D2.acc.diag
    var.D2.acc.diag
}
    
var.D2 <- function(t,k,eta,nA,nB)      # Does not include off-diagonal part of accordion!!!!!
{
    var.D2 <- var.D2.crabgrass(t,k,eta,nA,nB) + var.D2.acc.diag(t,k,eta,nA,nB)
    var.D2
}

# Computes the exact variance of D2
# Assumes periodic boundary conditions
#
# n: size of sequence A
# m: size of sequence B
# k: word size
# f: vector of words frequencies

varD2 <- function(n, m, k, f)
{
    p <- function(i, j=1)
    {
        sum(f^i)^j
    }
    f <- f / sum(f)

    sumVarYu <- n * m * (p(2, k) - p(2, 2 * k))

    # Non Uniform Case
    if (any(f != f[1]))
    {
        covCrab <- n * m * (n + m - 4 * k + 2) *
               (p(3, k) +
               2 * p(2, 2) * p(3) *
               (p(3, k - 1) - p(2, 2 * (k - 1))) /
               (p(3) - p(2, 2)) -
               (2 * k - 1) * p(2, 2 * k))
    }
    # Uniform Case
    else
    {
        covCrab <- 0
    }

    if (k == 1)
        return (sumVarYu + covCrab)

    covDiag  <- 2 * n * m *
            (p(2, k + 1) * (1 - p(2, k - 1)) / (1 - p(2)) - 
            (k - 1) * p(2, 2 * k))

    covAc1   <- 0
    for (t in 1:(k-1))
    {
        for (s in 0:(t - 1))
        {
            nu     <- floor((k - s) / (t - s))
            ro     <- (k -s) %% (t - s)
            covAc1 <- covAc1 +
                 (p(2, 2 * s) *
                 p(2 * nu + 3, ro) *
                 p(2 * nu + 1, t - s - ro) -
                 p(2, 2 * k))
        }
    }
    covAc1 <- covAc1 * 4 * n * m

    covAc2 <- 0
    r <- 1
    while (r < k)
    {
        t <- 1
        while (t < k)
        {
            nu    <- floor(k / (r + t))
            ro    <- k %% (r + t)
            prod1 <- 1
            for (i in 1:t)
            {
                l <- 1 + 2 * nu
                if (i <= ro)
                {
                    l <- l + 1
                }
                if (i <= (ro - r))
                {
                    l <- l + 1
                }
                prod1 <- prod1 * p(l)
            }
            prod2 <- 1
            for (i in 1:r)
            {
                l <- 1 + 2 * nu
                if (i <= ro)
                    l <- l + 1
                if (i <= (ro - t))
                    l <- l + 1
                prod2 <- prod2 * p(l)
            }
            covAc2 <- covAc2 + prod1 * prod2
            t <- t + 1
        }
        r <- r + 1
    }
    covAc2 <- covAc2 - ((k - 1)^2) * p(2, 2 * k)
    covAc2 <- covAc2 * 2 * n * m

    sumVarYu + covCrab + covDiag + covAc1 + covAc2
}

# Computes the exact mean of D2
# Assumes periodic boundary conditions
#
# n: size of sequence A
# m: size of sequence B
# k: word size
# f: vector of words frequencies

meanD2 <- function(n, m, k, f)
{
    f <- f / sum(f)
    n * m * sum(f^2)^k
}

main <- function()
{
    seq_sizes  <- c(10, 20, 40, 60, 100, 200, 400, 800)
    word_sizes <- c(1, 2, 4, 6, 8, 10)
    freqs      <- c(2, 2, 1, 1)
    baseDir    <- '/home/sf/reports/wordmatches/programs/results/distribution_fit/dna/non_uniform/periodic_results/'

    for (seq_size in seq_sizes)
    {
        for (word_size in word_sizes)
        {
            path <- paste(baseDir, 'D2_', seq_size, '_', word_size, sep='')
            d    <- NULL
            if (file.exists(path))
            {
                d <- scan(path, quiet=TRUE)
            }
            cat('n = m = ',   sprintf('%-4d', seq_size),
                ' ; k = ',    sprintf('%-3d', word_size),
                ' ; mean = ', sprintf('%.5e', meanD2(seq_size, seq_size, word_size, freqs)),
                ' ; var = ',  sprintf('%.5e', varD2(seq_size, seq_size, word_size, freqs)),
                #' ; sd = ',   sprintf('%.3e', var.D2(0, word_size, -1/3, seq_size, seq_size)),
                ifelse (is.null(d), '', paste(' ; est_mean =', sprintf('%.5e', mean(d)))),
                ifelse (is.null(d), '', paste(' ; est_var =',  sprintf('%.5e', var(d)))),
                '\n',
                sep='')
        }
    }
}

main()

# vim:ft=r:expandtab:ts=4:sw=4:sts=4:
