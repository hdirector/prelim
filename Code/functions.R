#Calculate periodogram
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

#Function to get autocorrelation for the matern
maternAc <- function(B, alpha, h, N, delta) {
  ac <- rep(NA, N)
  ac[1] <- ((B^2)*beta(0.5, alpha - 0.5))/(2*pi*abs(h)^(2*alpha  - 1)) #variance (cov at lag 0) 
  tau <- 1:(N - 1) #lags 
  ac[2:N] <- ((B^2*(abs(h)*delta*tau)^(alpha - 0.5)*besselK(abs(h)*delta*tau, alpha - 0.5))/
                (2^(alpha - 1/2)*pi^(1/2)*exp(lgamma(alpha))*abs(h)^(2*alpha  - 1))) 
  return(ac)
}

#Function to get autocorrelation for the OU process 
#(TO DO: DETERMINE IF TYPO IS IN THEIR CODE OR PAPER 
#(should there be a delta in this function)
ouAc <- function(A, w0, C, N, delta) {
  ac <- rep(NA, N)
  tau <- 0:(N - 1) #lags 
  ac <- (0i + A^2/abs(2*C))*exp(1i*w0*delta*tau)*exp(0i + -abs(C)*delta*tau)
  return(ac)
}

#Function to simulate a time series with complex Matern covariance
simMatern <- function(B, alpha, h, N, delta) {
  
  #calculate autocorrelation
  ac <- maternAc(B, alpha, h, N, delta)
  
  #Create covariance and generate realizations using standard gaussian tricks
  library(mvtnorm)
  L <- t(chol(toeplitz(ac)))
  u <- rmvnorm(N, mean = rep(0, 2), sigma = diag(2))
  Z <- L%*%u/sqrt(2) #renormalize variance
  colnames(Z) <- c("u", "v")
  return(Z)
}

#function to simulate a time series with fractional brownian motion
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

#Blurred Whittle likelihood
ll <- function(theta, Z, delta, curr, firstIndex, lastIndex, medIndex, incZero = FALSE, trans = TRUE) {
  #theta = (A, B, w0,c, h, alpha) 
  #A > 0: ou amplitude, B > 0: matern amplitude; w0: ou frequency, 
  #c > 0: ou dampening, h: matern slope, alpha: matern smoothness  (pg. 37) 
  if (trans == TRUE) {
    A <- expm1(theta[1]) + 1; B <- expm1(theta[2]) + 1; w0 <- theta[3]
    C <- expm1(theta[4]) + pi*sqrt(3)/N  + 1
    h <- expm1(theta[5]) + pi*sqrt(3)/N + 1
    alpha <- expm1(theta[6]) + 0.5 + 1
  } else {
    A <- theta[1]; B <- theta[2]; w0 <- theta[3];
    C <- theta[4]; h <- theta[5]; alpha <- theta[6]
  }
  
  #value is too high and will crash
  if (alpha > 171) {
    print("alpha")
  }
  N <- length(Z)
  tau <-  seq(0, N - 1, 1)
  sTau <- ouAc(A, w0, C, N, delta = 1) + maternAc(B, alpha, h, N, delta = 1)  
  sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) #Interpet this, why no negs for sBar?
  if (incZero == FALSE) {
    useIndex <- c(firstIndex:(medIndex - 1)) #check this line
    llVal <- sum(curr$sZ[useIndex]/sBar[useIndex] + log(sBar[useIndex]))
  } else {
    llVal <- sum(curr$sZ[firstIndex:lastIndex]/sBar[firstIndex:lastIndex] + log(sBar[firstIndex:lastIndex]))
  }
  print(llVal)
  return(100*log(llVal))
}

#Function to fit likelihood
fitModel <- function(Z, CF, delta, fracNeg, fracPos) {
  curr <- getPerio(Z, delta, dB = FALSE, noZero = FALSE)
  N <- length(Z)
  
  #set window of frequencies to consider in evalution (assuming frequencies sorted min to max)
  medIndex <- floor(N/2) + 1 #middle value (freq = 0)
  firstIndex <- round((medIndex - 1)*(1 - fracNeg) + 1) #index of minimum frequency considered
  lastIndex <- round(medIndex + fracPos*(N - medIndex)) #index of maximum frequency considered
  
  ####initial parameter estimates (using the well-reasoned parameters suggested in the paper's published code)
  #theta = (A, B, w0,c, h, alpha) 
  #A > 0: ou amplitude, B > 0: matern amplitude; w0: ou frequency, 
  #c > 0: ou dampening, h: matern slope, alpha: matern smoothness  (pg. 37) 
  parInit <- rep(NA, 6) 
  parInit[3] <- CF #set w0 to coriolis freq
  parInit[6] <- 1 #set smoothness to 1 (corresponds to f^-2 decay)
  #TO DO: UNDERSTAND LOGIC OF THESE INITIAL PARAMETERS (claim to be from solving simultaneous eq's of spectral, but not an exact match to what's in paper)
  
  #Find initial guess for OU amplitude and dampening, zero in on area around inertial oscillation peak 
  #Find index of  where frequencies move from inertial freq and turbulent background (as a guess, the author's use half the coriolis freq)
  divideIndex <- round(medIndex + (0.5*CF/pi)*medIndex)
  #Find index of max value (other than w = 0), need to search diff direction depending on sign of CF
  if (CF > 0) {
    maxIndex <- which.max(curr$sZ[divideIndex:N]) + divideIndex
  } else {
    maxIndex <- which.max(curr$sZ[1:divideIndex])
  }
  #Consider freq's one either side of peak for ou par (1/3 of the distance, arbitrary choice used by author's just for initial parameters)
  numTest <- ceiling(abs(medIndex - maxIndex)/3)
  ioIndex <- c((maxIndex - numTest):(maxIndex -1), (maxIndex + 1):(maxIndex + numTest))
  parInit[4] <- median((curr$sZ[ioIndex]*(delta*curr$omega[ioIndex] - CF)^2)/(curr$sZ[maxIndex] - curr$sZ[ioIndex]))
  parInit[1] <- curr$sZ[maxIndex]*parInit[4]
  parInit[1] <- sqrt(parInit[1]); parInit[4] <- sqrt(parInit[4])
  
  #Find initial guess for Matern amplitude and slope; zero in on area around turbulent background peak
  turbIndex <- c((medIndex - numTest):(medIndex - 1), (medIndex + 1):(medIndex + numTest))
  parInit[5] <- sqrt(median((curr$sZ[turbIndex]*(delta*curr$omega[turbIndex])^2)/(max(curr$sZ) - curr$sZ[turbIndex])))
  parInit[2] <- sqrt(max(curr$sZ))*parInit[5]
  
  #Transform parameters to an unconstrained space for optimization
  transParInit <- rep(NA, 6)
  transParInit[1] <- log1p(parInit[1] - 1) #0 < A < inf ===> -inf < log(A) < inf
  transParInit[2] <- log1p(parInit[2] - 1) #0 < B < inf ===> -inf < log(B) < inf 
  transParInit[3] <- parInit[3] #-inf < w0 < inf ===> -inf <- w0 < inf
  transParInit[4] <- log1p(parInit[4] - pi*sqrt(3)/N - 1) #pi*sqrt(3)/N < c < inf ===> -inf < log(c - pi*sqrt(3)/N) < inf
  transParInit[5] <- log1p(parInit[5] - pi*sqrt(3)/N - 1) #pi*sqrt(3)/N < c < inf ===> -inf < log(h - pi*sqrt(3)/N) < inf
  transParInit[6] <- log1p(parInit[6] - 0.5 - 1) #0.5 < alpha < inf ===> -inf < log(alpha - 0.5) < inf
  
  #Maximize likelihood numerically
  #theta = (A, B, w0,c, h, alpha) 
  #A > 0: ou amplitude, B > 0: matern amplitude; w0: ou frequency, 
  #c > 0: ou dampening, h: matern slope, alpha: matern smoothness  (pg. 37) 
  opt <- optim(transParInit, ll, Z = Z, delta = delta, curr = curr, 
               firstIndex = firstIndex, lastIndex = lastIndex, medIndex = medIndex,
               control = list(maxit = 1000000000, reltol=1e-1000))
  
  fin <-  c(expm1(opt$par[1]) + 1, expm1(opt$par[2]) + 1 , opt$par[3], expm1(opt$par[4]) + pi*sqrt(3)/N + 1, 
            expm1(opt$par[5]) + pi*sqrt(3)/N + 1, expm1(opt$par[6]) + 0.5 + 1)
  #fin <- truth
  return(list("A" = fin[1], "B" = fin[2], "w0" = fin[3],
              "C" = fin[4], "h" = fin[5], "alpha" = fin[6],
              "firstIndex" = firstIndex, "lastIndex" = lastIndex))

}