#Calculate periodogram
getPerio <- function(Z, delta, dB = TRUE) {
  library("waved")
  #fourier frquencies 
  N <- length(Z)
  if (N%%2 == 1) {
    omega <- (2*pi/(N*delta))*c(sort(-seq(1, (N - 1)/2, 1)), seq(1, (N - 1)/2, 1))
  } else {
    omega <- (2*pi/(N*delta))*c(sort(-seq(1, N/2, 1)), seq(1, N/2, 1))
  }
  #convert to decibals
  sZ <- (delta/N)*Mod(fft(Z))^2; sZ <- fftshift(sZ)
  if (dB == TRUE) {
    sZ <- 10*log10(sZ)
  }
  return(list("sZ" = sZ, "omega" = omega))
}

#Function to get autocorrelation for the matern
maternAc <- function(B, alpha, h, N) {
  ac <- rep(NA, N)
  ac[1] <- ((B^2)*beta(0.5, alpha - 0.5))/(2*pi*abs(h)^(2*alpha  - 1)) #variance (cov at lag 0) 
  tau <- 1:(N - 1) #lags 
  ac[2:N] <- ((B^2*(h*abs(tau))^(alpha - 0.5)*besselK(h*abs(tau), alpha - 0.5))/
                (2^(alpha - 1/2)*pi^(1/2)*gamma(alpha)*h^(2*alpha  - 1)))
  return(ac)
}

#Function to get autocorrelation for the OU process 
#(TO DO: DETERMINE IF TYPO IS IN THERE CODE OR PAPER 
#(should there be a delta in this function)
ouAc <- function(A, w0, c, tau, N) {
  ac <- rep(NA, N)
  tau <- 0:(N - 1) #lags 
  ac <- (0i + A^2/abs(2*c))*exp(1i*w0*tau)*exp(0i + -c*abs(tau))
  return(ac)
}

#Function to simulate a time series with complex Matern covariance
simMatern <- function(B, alpha, h, N) {
  
  #calculate autocorrelation
  ac <- maternAc(B, alpha, h, N)
  
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
ll <- function(theta, Z, delta) {
  A <- theta[1]; w0 <- theta[2]; c <- theta[3];
  B <- theta[4]; alpha <- theta[5]; h <- theta[6]
  N <- length(Z)
  tau <-  seq(0, N - 1, 1)
  curr <- getPerio(Z, delta)
  sTau <- ouAc(A, w0, c, tau, N) + maternAc(B, alpha, h, N)  
  sBar <- Re(fftshift(2*fft(sTau*(1 - tau/N)) - sTau[1])) 
  print(-sum(curr$sZ/sBar + log(sBar)))
  return(-sum(curr$sZ/sBar + log(sBar)))
}

