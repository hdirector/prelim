###Function to calculate periodogram###
#Inputs: Z: complex-valued velocity, dB: Boolean specifying whether result should be 
#returned in decibels or not, noZero: Boolean specifying whether to use the zero in 
#estimation or not 
#Outputs: list of omega (Fourier frequencies) and  sZ(spectral values at omega)
getPerio <- function(Z, delta, dB = TRUE, noZero = TRUE) {
  library("waved")
  #fourier frquencies 
  N <- length(Z)
  if(noZero == TRUE) {
    if (N%%2 == 1) {
      omega <- (2*pi/(N*delta))*c(sort(-seq(1, (N - 1)/2, 1)), seq(1, (N - 1)/2, 1))
    } else {
      omega <- (2*pi/(N*delta))*c(sort(-seq(1, N/2, 1)), seq(1, N/2, 1))
    }
  } else {
    if (N%%2 == 1) {
      omega <- (2*pi/(N*delta))*c(sort(-seq(1, (N - 1)/2, 1)), 0,  seq(1, (N - 1)/2, 1))
    } else {
      omega <- (2*pi/(N*delta))*c(sort(-seq(1, N/2, 1)), 0, seq(1, N/2, 1))
    }
  }
  #convert to decibals
  sZ <- (delta/N)*Mod(fft(Z))^2; sZ <- fftshift(sZ)
  if (dB == TRUE) {
    sZ <- 10*log10(sZ)
  }
  return(list("sZ" = sZ, "omega" = omega))
}

###Function to get autocorrelation for the Matern###
#Inputs: B: amplitude, alpha: smoothing, h: dampening, N: length of time series, delta: sampling interval
#Outputs: ac: autocovariance 
maternAc <- function(B, alpha, h, N, delta) {
  ac <- rep(NA, N)
  ac[1] <- ((B^2)*beta(0.5, alpha - 0.5))/(2*pi*abs(h)^(2*alpha  - 1)) #variance (cov at lag 0) 
  tau <- 1:(N - 1) #lags 
  ac[2:N] <- ((B^2*(abs(h)*tau)^(alpha - 0.5)*besselK(abs(h)*tau, alpha - 0.5))/
                (2^(alpha - 1/2)*pi^(1/2)*exp(lgamma(alpha))*abs(h)^(2*alpha  - 1))) 
  return(ac)
}

###Function to get autocorrelation for the OU process###
#Inputs: A: amplitude, w0: inertial frequency, C: dampening, N: length of time series, delta: sampling interval
#Outputs: ac: autocovariance 
ouAc <- function(A, w0, C, N, delta) {
  ac <- rep(NA, N)
  tau <- 0:(N - 1) #lags 
  ac <- (0i + A^2/abs(2*C))*exp(1i*w0*tau)*exp(0i + -abs(C)*tau)
  return(ac)
}

####Function to simulate a time series with Matern covariance###
#Inputs: B: amplitude, alpha: smoothing, h: dampening, N: length of time series, delta: sampling interval
#output: Z: sample time series with specified matern covariance
simMatern <- function(B, alpha, h, N, delta) {
  #calculate autocorrelation with Matern covariance
  ac <- maternAc(B, alpha, h, N, delta)
  #Create covariance and generate realizations using standard gaussian tricks
  library(mvtnorm)
  L <- t(chol(toeplitz(ac)))
  u <- rmvnorm(N, mean = rep(0, 2), sigma = diag(2))
  Z <- L%*%u/sqrt(2) #renormalize variance
  colnames(Z) <- c("u", "v")
  return(Z)
}

###Function to simulate a time series with fractional brownian motion(fbm)###
#Inputs: B: amplitude, alpha: smoothing, H: dampening, N: length of time series
#Outputs: Z: sample time series with specified fbm autocovariance 
simFBM <- function(B, alpha, H, N) {
  #Create covariance and generate realizations using standard gaussian tricks
  L <- matrix(nrow = N, ncol = N, data = NA)
  sigma2 <- (B^2)*(-gamma(2-2*(alpha - 0.5))*cos(pi*(alpha - 0.5))/(pi*(alpha - 0.5)*(2*alpha - 2)))
  H <- 2*alpha - 1
  for(s in 1:N) {
    for(t in 1:N) {
      #fraction brownian motion covariance
      L[t, s] <- (1/2)*(sigma2)*(t^H - abs(t -s)^H + s^H) 
    }
  }
  L <- t(chol(L))
  u <- rmvnorm(N, mean = rep(0, 2), sigma = diag(2))
  Z <- L%*%u/sqrt(2) #renormalize variance
  colnames(Z) <- c("u", "v")
  return(Z)
}

###Function for the Blurred Whittle likelihood with 6 parameters###
#Inputs: theta: vector of six parameters, theta = (A, B, w0,c, h, alpha) 
#A > 0: ou amplitude, B > 0: matern amplitude; w0: ou frequency, 
#c > 0: ou dampening, h > 0: matern slope, alpha > 1/2: matern smoothness  (pg. 37) 
#delta: sampling interval, curr: list of frequencies and spectral values (typically object 
#returned from getPerio function); firstIndex: first time point used; 
#lastIndex: last time point used; medIndex: time point of the 
#incZero: Boolean specifying if zero should be used in computation
#trans: Boolean indicating if optimization should be on the original constrained
#space or on a transformed unconstrained space
#output: negative of the log likelihood value 
ll <- function(theta, delta, curr, firstIndex, lastIndex, medIndex,
               incZero = FALSE, trans = TRUE) {
  N <- length(curr$sZ)
  #parameters to evaluate
  if (trans == TRUE) {
    A <- expm1(theta[1]) + 1; B <- expm1(theta[2]) + 1; w0 <- theta[3]
    C <- expm1(theta[4]) + pi*sqrt(3)/N  + 1
    h <- expm1(theta[5]) + pi*sqrt(3)/N + 1
    alpha <- expm1(theta[6]) + 0.5 + 1
  } else {
    A <- theta[1]; B <- theta[2]; w0 <- theta[3];
    C <- theta[4]; h <- theta[5]; alpha <- theta[6]
  }
  #calculate likelihood
  tau <-  seq(0, N - 1, 1)
  sTau <- ouAc(A, w0, C, N, delta = 1) + maternAc(B, alpha, h, N, delta = 1)  
  sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) #Interpet this, why no negs for sBar?
  if (incZero == TRUE) {
    llVal <- sum(curr$sZ[firstIndex:lastIndex]/sBar[firstIndex:lastIndex]) + sum(log(sBar[firstIndex:lastIndex]))
  } else {
    tempIndex <- firstIndex:lastIndex
    tempIndex <- tempIndex[-which(tempIndex == medIndex)]
    llVal <- sum(curr$sZ[tempIndex]/sBar[tempIndex]) + sum(log(sBar[tempIndex]))
  }
  return(llVal)
}

###Function for the Blurred Whittle likelihood with 5 parameters###
#Inputs: theta: vector of five parameters, theta = (A, B,c, h, alpha) 
#A > 0: ou amplitude, B > 0: matern amplitude; 
#c > 0: ou dampening, h > 0: matern slope, alpha > 1/2: matern smoothness  (pg. 37) 
#delta: sampling interval, curr: list of frequencies and spectral values (typically object 
#returned from getPerio function); firstIndex: first time point used; 
#lastIndex: last time point used; medIndex: time point of the 
#incZero: Boolean specifying if zero should be used in computation
#trans: Boolean indicating if optimization should be on the original constrained
#space or on a transformed unconstrained space
#output: negative of the log likelihood value 
llSimp <- function(theta, delta, curr, firstIndex, lastIndex, medIndex,
                   CF, incZero = FALSE, trans = TRUE) {
  N <- length(curr$sZ)
  #parameters to evaluate
  if (trans == TRUE) {
    A <- expm1(theta[1]) + 1; B <- expm1(theta[2]) + 1; 
    C <- expm1(theta[3]) + pi*sqrt(3)/N  + 1
    h <- expm1(theta[4]) + pi*sqrt(3)/N + 1
    alpha <- expm1(theta[5]) + 0.5 + 1
  } else {
    A <- theta[1]; B <- theta[2]; 
    C <- theta[3]; h <- theta[4]; alpha <- theta[5]
  }
  #calculate likelihood
  tau <-  seq(0, N - 1, 1)
  sTau <- ouAc(A, CF, C, N, delta = 1) + maternAc(B, alpha, h, N, delta = 1)  
  sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
  if (incZero == TRUE) {
    llVal <- sum(curr$sZ[firstIndex:lastIndex]/sBar[firstIndex:lastIndex]) + sum(log(sBar[firstIndex:lastIndex]))
  } else {
    tempIndex <- firstIndex:lastIndex
    tempIndex <- tempIndex[-which(tempIndex == medIndex)]
    llVal <- sum(curr$sZ[tempIndex]/sBar[tempIndex]) + sum(log(sBar[tempIndex]))
  }
  return(llVal)
}


###Function to fit likelihood###
#Inputs: Z: complex-value velocity vector, CF: coriolis frequency
#delta: sampling increment
#fracNeg: percent of the negative values used in model fitting
#fracPos: percent of the positive value used in model fitting
#quantSet: quantile used when picking initial parameters
#incZero: Boolean indicating where 
#hess: default value for Hessian
#needInits: Boolean indicating whether initial values need to be 
#calculated for optimization or not
#parInit: initial values that can be supplied
#getHess: Boolean indicating whether Hessian needs to be calculated or not
#simpModel: Boolean indicating that the simplified 5 parameter model should be fit rather 
#than the full six parameter model
#Outputs: List of parameter estimates (5 or 6 values depending on whether fitting full or 
#simplified model, firstIndex: first time point used in fitting, lastIndex: last time point 
#used in fitting; hess: hessian; llVal: negative of the log likelihood value
fitModel <- function(Z, CF, delta, fracNeg, fracPos, quantSet, incZero = FALSE, hess = NULL,
                     needInits = TRUE, parInit = NULL, getHess = FALSE, simpModel = FALSE) {
  
  #calculate periodogram for data
  curr <- getPerio(Z, delta, dB = FALSE, noZero = FALSE)
  
  N <- length(curr$sZ)
  
  #set window of frequencies to consider in evalution (assuming frequencies sorted min to max)
  medIndex <- floor(N/2) + 1 #middle value (freq = 0)
  firstIndex <- round((medIndex - 1)*(1 - fracNeg) + 1) #index of minimum frequency considered
  lastIndex <- round(medIndex + fracPos*(N - medIndex)) #index of maximum frequency considered
  
  if (needInits == TRUE) {
    parInit <- rep(NA, 6) 
    parInit[3] <- CF #set w0 to coriolis freq
    parInit[6] <- 1 #set smoothness to 1 (corresponds to f^-2 decay)
    #TO DO: UNDERSTAND LOGIC OF THESE INITIAL PARAMETERS (claim to be from solving simultaneous eq's of spectral, but not an exact match to what's in paper)
    
    #Find initial guess for OU amplitude and dampening, zeroing in on area around inertial oscillation peak 
    #Find index of  where frequencies move from inertial freq and turbulent background (as a guess, Sykulski et al. use half the coriolis freq)
    divideIndex <- round(medIndex + (0.5*CF/pi)*medIndex)
    #Find index of max value (other than frequency = 0); need to search different direction depending on 
    #sign of CF (positive in Southern Hemisphere and negative in Northern Hemisphere)
    if (CF > 0) {
      maxIndex <- which.max(curr$sZ[divideIndex:N]) + divideIndex
    } else {
      maxIndex <- which.max(curr$sZ[1:divideIndex])
    }
    #Consider freq's one either side of peak for OU parameter (Use 1/3 of the distance following Sykulski et al. )
    numTest <- floor(abs(medIndex - maxIndex)/3)
    ioIndex <- c((maxIndex - numTest):(maxIndex -1), (maxIndex + 1):(maxIndex + numTest))
    parInit[4] <- quantile((curr$sZ[ioIndex]*(delta*curr$omega[ioIndex] - CF)^2)/(curr$sZ[maxIndex] - curr$sZ[ioIndex]), quantSet)
    parInit[1] <- curr$sZ[maxIndex]*parInit[4]
    parInit[1] <- sqrt(parInit[1]); parInit[4] <- sqrt(parInit[4])
    
    #Find initial guess for Matern amplitude and slope. Zero in on area around turbulent background peak
    turbIndex <- c((medIndex - numTest):(medIndex - 1), (medIndex + 1):(medIndex + numTest))
    parInit[5] <- quantile(((curr$sZ[turbIndex]*(delta*curr$omega[turbIndex])^2)/(max(curr$sZ) - curr$sZ[turbIndex])), quantSet)
    parInit[2] <- sqrt(max(curr$sZ)*parInit[5]); parInit[5] <- sqrt(parInit[5])
  }  
  
  #Transform parameters to an unconstrained space for optimization
  if (simpModel == FALSE) {
    transParInit <- rep(NA, 6)
    transParInit[1] <- log1p(parInit[1] - 1) #0 < A < inf ===> -inf < log(A) < inf
    transParInit[2] <- log1p(parInit[2] - 1) #0 < B < inf ===> -inf < log(B) < inf 
    transParInit[3] <- parInit[3] #-inf < w0 < inf ===> -inf <- w0 < inf
    transParInit[4] <- log1p(parInit[4] - pi*sqrt(3)/N - 1) #pi*sqrt(3)/N < c < inf ===> -inf < log(c - pi*sqrt(3)/N) < inf
    transParInit[5] <- log1p(parInit[5] - pi*sqrt(3)/N - 1) #pi*sqrt(3)/N < h < inf ===> -inf < log(h - pi*sqrt(3)/N) < inf
    transParInit[6] <- log1p(parInit[6] - 0.5 - 1) #0.5 < alpha < inf ===> -inf < log(alpha - 0.5) < inf
  } else {
    if (needInits == TRUE) {
      parInit <- parInit[c(1:2, 4:6)]
    }
    transParInit <- rep(NA, 5)
    transParInit[1] <- log1p(parInit[1] - 1) #0 < A < inf ===> -inf < log(A) < inf
    transParInit[2] <- log1p(parInit[2] - 1) #0 < B < inf ===> -inf < log(B) < inf 
    transParInit[3] <- log1p(parInit[3] - pi*sqrt(3)/N - 1) #pi*sqrt(3)/N < c < inf ===> -inf < log(c - pi*sqrt(3)/N) < inf
    transParInit[4] <- log1p(parInit[4] - pi*sqrt(3)/N - 1) #pi*sqrt(3)/N < h < inf ===> -inf < log(h - pi*sqrt(3)/N) < inf
    transParInit[5] <- log1p(parInit[5] - 0.5 - 1) #0.5 < alpha < inf ===> -inf < log(alpha - 0.5) < inf
  }
  
  #Maximize likelihood numerically
  library("pracma")
  if (simpModel == FALSE) {
    opt <- fminsearch(ll, transParInit, delta = delta, curr = curr,firstIndex = firstIndex,
                       lastIndex = lastIndex, medIndex = medIndex, maxiter = 3000) 
  } else {
    opt <- fminsearch(llSimp, transParInit, delta = delta, curr = curr, 
                      firstIndex = firstIndex, lastIndex = lastIndex, CF = CF,
                      medIndex = medIndex, maxiter = 3000)
  }
  llVal <- opt$fval
  
  #Convert parameters back to natural scale
  if (simpModel == FALSE) {
    fin <-  c(expm1(opt$xval[1]) + 1, expm1(opt$xval[2]) + 1 , opt$xval[3], expm1(opt$xval[4]) + pi*sqrt(3)/N + 1, 
              expm1(opt$xval[5]) + pi*sqrt(3)/N + 1, expm1(opt$xval[6]) + 0.5 + 1)
  } else {
    fin <-  c(expm1(opt$xval[1]) + 1, expm1(opt$xval[2]) + 1, expm1(opt$xval[3]) + pi*sqrt(3)/N + 1, 
              expm1(opt$xval[4]) + pi*sqrt(3)/N + 1, expm1(opt$xval[5]) + 0.5 + 1)
  }
  
  #Calculate Hessian where needed
  if (simpModel == FALSE & getHess == TRUE) {
    library("optimx")
    hess <- optimHess(fin, ll,  delta = delta, curr = curr, trans = FALSE,
                  firstIndex = firstIndex, lastIndex = lastIndex, medIndex = medIndex,
                  control = list(maxit = 2000))
  } 
  
  #Return results
  if (simpModel == FALSE) {
    return(list("A" = fin[1], "B" = fin[2], "w0" = fin[3],
                "C" = fin[4], "h" = fin[5], "alpha" = fin[6],
                "firstIndex" = firstIndex, "lastIndex" = lastIndex,
                "hess" = hess, llVal = llVal))
  } else {
    return(list("A" = fin[1], "B" = fin[2], 
                "C" = fin[3], "h" = fin[4], "alpha" = fin[5],
                "firstIndex" = firstIndex, "lastIndex" = lastIndex,
                "hess" = hess, llVal = llVal))
  }
  
}