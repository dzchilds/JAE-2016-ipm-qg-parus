
library(inline)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to handle the numerical implementation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate mesh delta -> 'x' 
#
mkdx <- function(iPar) {
  (iPar$maxx-iPar$minx)/iPar$msizex
}

# calculate mesh delta -> 'g' 
#
mkdg <- function(iPar) {
  (iPar$maxg-iPar$ming)/iPar$msizeg
}

# construct the integration mesh -> 'x' 
# 
mkmeshx <- function(iPar) {
  iPar$minx+mkdx(iPar)*(1:iPar$msizex-1/2)
}

# construct the integration mesh -> 'g' 
#
mkmeshg <- function(iPar) {
  iPar$ming+mkdg(iPar)*(1:iPar$msizeg-1/2)
}

# calculate the control parameters
# 
calcControlParams <- function(iPar)
{
  # extract the max age
  maxA <- iPar$maxA
  if (maxA<3) stop("Must track at least 3 age classes")
  # derive the 'delta'
  dx <- mkdx(iPar); dg <- mkdg(iPar)
  # derive the mesh points for quadrature (!! and remember the centering !!)
  meshx <- mkmeshx(iPar); meshg <- mkmeshg(iPar)
  # number of mesh points
  nx <- iPar$msizex
  ng <- iPar$msizeg
  # construct the 'grid' values for the state variable approximations
  grid.nt <- list(F=NULL, M=NULL)
  grid.nt $ F <- as.list(expand.grid(x=meshx, g=meshg, A=seq.int(1,maxA)))
  attributes(grid.nt $ F) <- NULL
  names(grid.nt $ F) <- c("x","g","A")
  grid.nt $ M <- as.list(expand.grid(g=meshg, A=seq.int(1,maxA)))
  attributes(grid.nt $ M) <- NULL
  names(grid.nt $ M) <- c("g","A")
  # construct the age indexing variable for the state variables
  iA.nt <- list(F=list(), M=list())
  for (A in seq.int(1,maxA)) {
    iA.nt $ F [[A]] <- which(grid.nt $ F $ A == A)
    iA.nt $ M [[A]] <- which(grid.nt $ M $ A == A)
  }
  # return all this
  cPar <- list(nx = nx, dx = dx, meshx = meshx,
               ng = ng, dg = dg, meshg = meshg,
               iA.nt = iA.nt, grid.nt = grid.nt)
  return(c(iPar, cPar))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions defining the vital rates 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract the correct intercept for the synchrony model 
# -- shared environmental component of the trait
#
evalIntercept <- function(mPar, A)
{
  as.vector(mPar[c("pLD.iAf0", "pLD.iAf1")])[ifelse(A == 1, 1, 2)]
}

# survival function
# -- parameterised in terms of the 'observed phenotype'
#    i.e. breeding value + residual + intercept = g + x + I
#
f.SU.F <- function(grid.nt, Nt, Yr, Pe, mPar)
{
    # construct centred / derived variables
    xg <- grid.nt $ x  + grid.nt $ g + evalIntercept(mPar, grid.nt $ A)          
    A  <- grid.nt $ A  - mPar["cenA"]
    Nt <- Nt - mPar["cenN"]
    Yr <- Yr - mPar["cenY"]
    # compute and return expected survival
    pNames <-
        paste("pSU",
              c("(Intercept)","layDt","I(layDt^2)","Ag","I(Ag^2)","Nt","Yr","Pe"),
              sep=".")
    p  <- mPar[pNames]
    nu <- p[1] + p[2]*xg + p[3]*xg^2 + p[4]*A + p[5]*A^2 + p[6]*Nt + p[7]*Yr +p[8]*Pe
    return(1/(1+exp(-nu)))
}

# recruitment function ('FC' in name references fecundity)
# -- parameterised in terms of the 'observed phenotype'
#    i.e. breeding value + residual + intercept = g + x + I
#
f.FC.F <- function(grid.nt, Nt, Yr, Pe, mPar)
{
    # construct centred / derived variables
    xg <- grid.nt $ x  + grid.nt $ g + evalIntercept(mPar,  grid.nt $ A)
    A  <- grid.nt $ A  - mPar["cenA"]
    Nt <- Nt - mPar["cenN"]
    Yr <- Yr - mPar["cenY"]
    # compute and return expected recruits
    pNames <-
        paste("pFC",
              c("(Intercept)","layDt", "I(layDt^2)","Ag","I(Ag^2)","Nt","Yr","Pe"),
              sep=".")
    p  <- mPar[pNames]
    nu <- p[1] + p[2]*xg + p[3]*xg^2 + p[4]*A + p[5]*A^2 + p[6]*Nt + p[7]*Yr +p[8]*Pe 
    return(exp(nu))
}

# 'mother-daughter' covariance function
# -- associated with the autocor ** residual ** component of variation (i.e. the l-state)
# -- no maternal effect here
#
f.MD.F <- function(x_, mPar)
{
    pNames <- paste("pLD", c("iAf0", "sigmaR"), sep=".")
    p <- mPar[pNames]
    return(dnorm(x_, 0, p[2]))
}

# 'ontogeny' function
# -- associated with the autocor ** residual ** component of variation (i.e. the l-state)
#
f.ON.F <- function(x_, x, mPar)
{
    pNames <- paste("pLD", c("iAf1","sigmaR", "corr"), sep=".")
    p <- mPar[pNames]
    mu <- p[3]*x
    sg <- sqrt( (1-p[3]^2)*p[2]^2 )
    return(dnorm(x_, mu, sg))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convenience functions to evaluate marginal distributions of the state distribution 
# -- use 'with' inside a function (naughty, but safe here)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# evaluate marginal distribution -> females
#
getMargF <- function(f, cPar, type)
with(cPar, {
    dim(f) <- c(nx, ng, maxA)
    dimnames(f) <- list(x=NULL, g=NULL, A=NULL)
    fmarg <- apply(f, type, sum)
    deltas <- c(dx, dg, 1)
    fmarg <- fmarg * prod(deltas) / prod(deltas[match(type, c("x","g","A"))], na.rm=TRUE)
    return(fmarg)
})

# evaluate marginal distribution -> males
# 
getMargM <- function(f, cPar, type)
with(cPar, {
    dim(f) <- c(ng, maxA)
    dimnames(f) <- list(g=NULL, A=NULL)
    fmarg <- apply(f, type, sum)
    deltas <- c(dg, 1)
    fmarg <- fmarg * prod(deltas) / prod(deltas[match(type, c("g","A"))], na.rm=TRUE)
    return(fmarg)
})

# wrapper so we only have to call one function
# 
getMarg <- function(sx, f, cPar, type) {
    if (sx=="F") return(getMargF(f, cPar, type))
    if (sx=="M") return(getMargM(f, cPar, type))
    stop ("invalid sex indicator")
}

# C++ code for 1d convolution
# -- need a configured compiler for this to work (fine on osx / linux, maybe windows too?)
#
convolve_1d <- cxxfunction(signature(d1S = "numeric", d2S = "numeric", dxS = "numeric"),
  plugin = "Rcpp", '
    Rcpp::NumericVector d1(d1S), d2(d2S);
    int x_d1 = d1.length(), x_d2 = d2.length();
    double dx = Rcpp::as<double>(dxS);

    Rcpp::NumericVector output(x_d1 + x_d2 - 1);
    for (int p = 0; p < x_d1; p++) {
        for (int i = 0; i < x_d2; i++) {
            output(p + i) += d1(p) * d2(i) * dx;
        }
    }
    return output;
')

# calculate the offspring breedng value distribution from dam / sire distributions  
# -- only designed to work if the breeding value mesh is ** centered at zero **
# -- read the previous line!
#
getOffgDistEvo <- function(phiF, phiM, meshg, Va)
{
    # makes we use the mesh to define min/max etc
    ming <- min(meshg)
    maxg <- max(meshg)
    ng <- length(meshg)
    dg <- diff(meshg)[1]
    # convolution female g/2 dist with male g/2 dist
    meshg2 <- seq(ming, maxg, length=2*ng-1)
    phi <- 2 * convolve_1d(phiF, phiM, dg)
    dg2 <- diff(meshg2)[1]
    # convolution of resultant dist with Va/2 normal
    meshg3 <- seq(2*ming, 2*maxg, length=4*ng-3)
    phi <- convolve_1d(phi, dnorm(meshg2, sd=sqrt(Va/2)), dg2)
    # return the submesh corresponding to the original mesh
    ii <- seq.int(ng, 3*ng-1, by=2)
    return(phi[ii])
}

# calculate the offspring breedng value distribution using mother-offspring regression
# -- use this to implement Yngvild & Langangen (2015)
# -- not used here or in paper, just for interest really
# -- doesn't evolve, implementation not correct?..
#
getOffgDistReg <- function(phiF, phiM, meshg, Va)
{
  return(dnorm(meshg, mean = phiF / 2, sd = sqrt(3*Va/4)))
}

# calculate the offspring breedng value distribution assuming no evolution
# -- i.e., 'g' is just another non-heritable source of individual variation
# 
getOffgDistNoE <- function(phiF, phiM, meshg, Va)
{
    return(dnorm(meshg, mean=0, sd=sqrt(Va)))
}

# construct initial state distribution
# -- this is a list, not an atomic vector
#
makeInitState <- function(cPar, mPar, Nt, sigmaG)
{
  with(cPar, {
    ii <- unlist(iA.nt $ F [1])
    ntF <- numeric(nx * ng * maxA)
    ntF[ii] <- Nt * dnorm(grid.nt$F$x[ii], mean = 0, sd = 7) *
      dnorm(grid.nt$F$g[ii], mean = 0, sd = 2)
    ii <- unlist(iA.nt $ M [1])
    ntM <- numeric(ng * maxA)
    ntM[ii] <- Nt * dnorm(grid.nt$M$g[ii], mean = 0, sd = 2)
    return(list(F=ntF, M=ntM))
  })
}

#  iterate the system
#
iterateEcoEvoDynamics <- function(mPar, cPar, eSet, nt0, keepAll=FALSE) with(c(cPar,eSet),
{
    nSim <- length(ySet)
    # make sure we have vectors where they are expected
    if (length(yTrend)==1) {
        yTrend <- rep(yTrend, nSim)
        eSet$yTrend <- yTrend
    }
    if (length(period)==1) {
        period <- rep(period, nSim)
        eSet$period <- period
    }
    # initial conditions
    ntF <- nt0 $ F; ntM <- nt0 $ M
    # storage lists
    ntFts <- list(ntF); ntMts <- list(ntM)
    # now run the simulation
    for (tt in seq.int(1, nSim)) {

        # extract the current year type
        mParNow <- mPar[, ySet[tt]]
        # extract temporal trend terms
        Yr <- yTrend[tt]; Pe <- period[tt]

        # 1) *** fecundity ***
        # current density of breeding females 
        NtF <- sum(ntF * dx * dg)
        # post-fertility density function (distribution of maternal traits in recruits)
        fn <- f.FC.F(grid.nt$F, NtF, Yr, Pe, mParNow) * ntF
        # total number of new recruits
        Nt1 <- sum(fn) * dx * dg
        # normalised post-fertility density function
        delF <- fn / Nt1
        # marginal maternal genotype distribution
        phiF <- getMargF(delF, cPar, type="g")
        # normalised male post-fertility density function
        delM <- ntM / (sum(ntM) * dg) 
        # marginal paternal genotype distribution
        phiM <- getMargM(delM, cPar, type="g")
        # get the offspring genotype distribution
        omega <- getOffgDist(phiF, phiM, meshg, Va=mParNow["pLD.sigmaG"]^2)
        # construct the joint distribution of 'x' and 'g' in female recruits
        ntF1 <- outer(meshx, omega, function(x_, omega) omega * f.MD.F(x_, mParNow))
        dim(ntF1) <- NULL
        ntF1 <- Nt1 * ntF1 / 2
        # construct the distribution 'g' in male recruits
        ntM1 <- omega
        ntM1 <- Nt1 * ntM1 / 2
        
        # 2) *** survival ***
        # marginal age distribution, pre-survival
        preNtFbyA <- getMargF(ntF, cPar, type="A")
        # post-survival density function of females
        fn <- f.SU.F(grid.nt$F, NtF, Yr, Pe, mParNow) * ntF
        # 'ontogeny'
        gn <- outer(meshx, meshx, f.ON.F, mParNow)
        dim(fn) <- c(nx, ng * maxA)
        fn <- gn %*% fn * dx
        dim(fn) <- NULL
        # marginal age distribution, post-survival
        posNtFbyA <- getMargF(fn, cPar, type="A")
        # age-specific survival probability
        pSurv <- sum(posNtFbyA) / sum(preNtFbyA)
        # pSurv <- ifelse(is.nan(pSurv), 0, pSurv)

        # 3) *** move everyone around among age classes *** 
        # survivors
        for (A in seq.int(maxA-1, 1)) {
            ntF[ unlist(iA.nt$F[A+1]) ] <-  fn[ unlist(iA.nt$F[A]) ]
            ntM[ unlist(iA.nt$M[A+1]) ] <- ntM[ unlist(iA.nt$M[A]) ] * pSurv
        }
        # new recruits
        ntF[ unlist(iA.nt$F[1]) ] <- ntF1
        ntM[ unlist(iA.nt$M[1]) ] <- ntM1

        # 4) *** store everything ***
        ntFts[[ tt+1 ]] <- ntF
        ntMts[[ tt+1 ]] <- ntM
    }
    return(list(F=ntFts, M=ntMts, eSet=eSet, mPar=mPar))
})


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to compute terms of the age-structured Price (ASPE) equation 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# evaluate expectation of some function ('f') wrt to state distribution -> females
# -- 'f' must be an expression involving 'g' and/or 'x'
#
calcExpectationF <- function(n.A.z, f, mPar, cPar, byA=TRUE) {
    with(cPar, {
        wgts <- dx * dg
        fx <- eval(f, envir = grid.nt$F)
        integrands <-  fx * n.A.z
        if (byA==TRUE) {
            dim(integrands) <- c(nx, ng, maxA)
            return(apply(integrands, 3, sum) * wgts / getMargF(n.A.z, cPar, type="A"))
        } else {
            return(sum(integrands) * wgts)
        }
    })
}

# evaluate expectation of some function ('f') wrt to state distribution -> males
# -- 'f' must be an expression involving 'g' and/or 'x'
#
calcExpectationM <- function(n.A.z, f, mPar, cPar, byA=TRUE) {
    with(cPar, {
        wgts <- dg
        fx <- eval(f, envir = grid.nt$M)
        integrands <-  fx * n.A.z
        if (byA==TRUE) {
            dim(integrands) <- c(ng, maxA)
            return(apply(integrands, 2, sum) * wgts / getMargM(n.A.z, cPar, type="A"))
        } else {
            return(sum(integrands) * wgts)
        }
    })
}

# convenience function to call the sex-specific expectation function above
#
calcExpectation <- function(sx, n.A.z, f, mPar, cPar, byA=TRUE) {
    if (sx=="F") return(calcExpectationF(n.A.z, f, mPar, cPar, byA=byA))
    if (sx=="M") return(calcExpectationM(n.A.z, f, mPar, cPar, byA=byA))
    stop ("Invalid sex indicator")
}

# calculate the various parts (means, weights, etc) of ASPE 
# -- works with a long time series produced by 'iterateEcoEvoDynamics'
# 
calcPriceDecompTerms <- function(cPar, simSys, vrType="g") with(c(cPar, simSys),
{

    ySet   <- eSet $ ySet
    yTrend <- eSet $ yTrend
    period <- eSet $ period
    
    nSim <- length(ySet) 
    priceTerms <- list()
    # set the expectation we will work with for females
    expectF <- switch(vrType,
                      g = expression( g ),
                      p = expression( x + g + evalIntercept(mPar, A) ),
                      stop("Did not recognise variable"))
    expectM <- switch(vrType,
                      g = expression( g ),
                      p = expression( return(NA) ),
                      stop("Did not recognise variable"))
        
    # iterate over the time series
    for (tt in seq.int(1, nSim)) {

        yrType <- ySet[tt]
        # extract the current year type
        mParNow <- mPar[, yrType]
        # extract temporal trend terms
        Yr <- yTrend[tt]       
        Pe <- period[tt]
        #
        N <- N.p1 <- c(F = 0, M = 0)
        phi <- n3.A.g <- list(F=NULL, M=NULL)
        priceTerms [[tt]] <- list(F = NULL, M = NULL, T = NULL)
        for (sx in c("M", "F")) {
            
            # population vectors
            n.A.z    <- simSys [[ sx ]] [[ tt   ]] 
            n.A.z.p1 <- simSys [[ sx ]] [[ tt+1 ]]
            # number in each age class this year
            N.A    <- getMarg(sx, n.A.z,    cPar, type="A") 
            # number in each age class next year
            N.A.p1 <- getMarg(sx, n.A.z.p1, cPar, type="A")
            # density now, density next year
            N[sx]    <- sum(N.A)
            N.p1[sx] <- sum(N.A.p1)
            if (sx=="F") {
                # mean conditional on age
                Z.1 <- calcExpectation(sx, n.A.z,  expectF, mParNow, cPar)
                # survival function
                survFunc <- f.SU.F(grid.nt$F, N[sx], Yr, Pe, mParNow)
                # age-specific survival probs
                S.A <- calcExpectation(sx, n.A.z, survFunc, mParNow, cPar)
                # post survival mean conditional on age
                n2.A.z <- survFunc * n.A.z
                Z.2 <- calcExpectation(sx, n2.A.z, expectF, mParNow, cPar)
                # recruitment function
                recrFunc <- f.FC.F(grid.nt$F, N[sx], Yr, Pe, mParNow)
                # age-specific recruitment rates ! don't forget to scale by sex ratio !
                R.A <- calcExpectation(sx, n.A.z, recrFunc, mParNow, cPar) / 2
                # post fecundity mean conditional on age
                n3.A.z <- recrFunc * n.A.z
                Z.3 <- calcExpectation(sx, n3.A.z, expectF, mParNow, cPar)
                # post ontogeny mean conditional on age !!! HERE IS YOUR MISTAKE - forgot offset !!!
                Z.4 <- calcExpectation(sx, n.A.z.p1, expectF, mParNow, cPar)
                Z.4 <- c(Z.4[-1], NA)
                if (vrType=="p")
                    Z.4 <- Z.4 + mPar[ "pLD.iAf1", ySet[tt+1] ] -  mPar[ "pLD.iAf1", ySet[tt] ]
            } else if (sx=="M") {
                Z.1 <- calcExpectation(sx, n.A.z, expectM, mParNow, cPar)
                Z.4 <- Z.3 <- Z.2 <- Z.1
                Z.4[length(Z.4)] <- NA  # not really needed but keeps F and M data the same
                S.A <- c(N.A.p1[-1] / N.A[-length(N.A)], 0)
                R.A <- rep( N.A.p1[1] / sum( N.A ), maxA)
                n3.A.z <- n.A.z
            } else stop("Unrecognised sex indicator/n")
            # marginal "g x A" and "g" distribution of breeding individual's genotypes
            C.norm <- getMarg(sx, n3.A.z, cPar, type="A")
            n3.A.g [[sx]] <- getMarg(sx, n3.A.z, cPar, type=c("A","g")) / C.norm
            phi [[sx]] <- getMarg(sx, n3.A.z, cPar, type="g") / sum(C.norm)    
            # store the important quantities ! this generates warnings because R wants to grab
            # rownames from the subsetted 'N' and 'N.p1' vectors which are length 1
            priceTerms [[tt]] [[sx]] <-
                suppressWarnings(data.frame(N = N[sx], N.A = N.A, N.p1 = N.p1[sx],
                                            W = N.p1[sx] / N[sx], S.A = S.A, R.A = R.A,
                                            c.A = N.A / N[sx], c.A.p1 = c(N.A.p1[-1], 0) / N.p1[sx],
                                            Z.1 = Z.1, Z.2 = Z.2, Z.3 = Z.3, Z.4 = Z.4, Z.5 = NA))
        }
        #  get the offspring genotype distribution
        Va <- mParNow["pLD.sigmaG"] ^ 2
        # offspring genotype distribution of male/female parent, by age
        omegaF <- apply(n3.A.g[["F"]], 1, function(phiF) getOffgDist(phiF, phi[["M"]], meshg, Va))
        omegaM <- apply(n3.A.g[["M"]], 1, function(phiM) getOffgDist(phi[["F"]], phiM, meshg, Va))
        # obtain offspring genoptype means
        Z.5.F <- apply(omegaF, 2, function(f) sum(meshg * f) * dg)
        Z.5.M <- apply(omegaM, 2, function(f) sum(meshg * f) * dg)
        # add the environmental component if required
        if (vrType=="p") {
            Z.5.F <- Z.5.F + mPar[ "pLD.iAf0", ySet[tt+1] ]
            Z.5.M <- rep(NA, maxA)         
        }
        # store the offspring trait means
        priceTerms [[tt]] $ F [,"Z.5"] <- Z.5.F
        priceTerms [[tt]] $ M [,"Z.5"] <- Z.5.M
        # add the total population metrics
        priceTerms [[tt]] $ T <- list(q    = N["F"] / sum(N),
                                      q.p1 = N.p1["F"] / sum(N.p1),
                                      W    = sum(N.p1) / sum(N))
    }
    return(priceTerms)
})

# calculate the main terms of the one-sex ASPE 
# -- works with a long time series produced by 'calcPriceDecompTerms'
# 
calcPrice1Sx <- function(priceTerms, nfig=16)
    within(priceTerms, {
        nA <- length(c.A)
        SD.S <- round( (c.A * S.A / W) * (Z.2 - Z.1), nfig )
        SD.F <- round( (c.A * R.A / W) * (Z.3 - Z.1), nfig )
        ON   <- round( (c.A * S.A / W) * (Z.4 - Z.2), nfig )
        PO   <- round( (c.A * R.A / W) * (Z.5 - Z.3), nfig )
        AF.F <- round( (c.A * R.A / W) * Z.1, nfig)
        AF.1 <- round( (c.A.p1 * Z.1), nfig)
        AF.0 <- round( - c.A * Z.1, nfig)
        SD.S[nA] <- 0
        ON  [nA] <- 0
        AF.1[nA] <- 0
        rm(nA)
    }
)

# calculate the main terms of the two-sex ASPE
# -- works with a long time series produced by 'calcPriceDecompTerms'
# 
calcPrice2Sx <- function(priceTerms, nfig=128)
{
    priceTerms <- within(priceTerms, {
        with(T, {
            envp <- parent.env(environment())
            assign("scF", q / W,     pos=envp)
            assign("qF" , q,         pos=envp)
            assign("scM", (1-q) / W, pos=envp)
            assign("qM" , (1-q),     pos=envp)
        })
        F <- within(calcPrice1Sx(F, nfig=nfig), {
            SD.S <- SD.S * scF * W; SD.F <- SD.F * scF * W
            ON   <- ON   * scF * W; PO   <- PO   * scF * W
            AF.1 <- AF.1 * scF * W; AF.0 <- AF.0 * qF; AF.F <- AF.F * scF * W
        })
        M <- within(calcPrice1Sx(M, nfig=nfig), {
            SD.S <- SD.S * scM * W; SD.F <- SD.F * scM * W
            ON   <- ON   * scM * W; PO   <- PO   * scM * W
            AF.1 <- AF.1 * scM * W; AF.0 <- AF.0 * qM; AF.F <- AF.F * scM * W
        })  
        rm(scF, qF, scM, qM)
    })
    return(priceTerms)
}