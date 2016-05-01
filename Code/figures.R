###Hannah Director###
#Research paper prelim 
rm(list = ls())
source("/users/hdirector/Documents/prelim/prelim/Code/functions.R")

#Constants
DELTA <- 2

#Convert data from matlab to R
library("R.matlab")
drifterbetty2 <- readMat("/users/hdirector/Documents/prelim/prelim/Code/drifterbetty2.mat")
drift <- list("num" = as.vector(drifterbetty2$drifterbetty[,,1]$num),
              "lat" = as.vector(drifterbetty2$drifterbetty[,,1]$lat),
              "lon" = as.vector(drifterbetty2$drifterbetty[,,1]$lon),
              "cv" = as.vector(drifterbetty2$drifterbetty[,,1]$cv),
              "f" = as.vector(drifterbetty2$drifterbetty[,,1]$f))


###Figure 1 (right)###
#Indices of sections of interest
startTs <- 1901; endTs <- 4900 #index of period of interest is obtained from Sykulski et al posted Matlab code
bStartTs <- 2971; bEndTs <-3570 #index of 50 day period in blue is obtained from Sykulski et al posted Matlab code
gStartTs <- 4201; gEndTs <- 4800 #index of 50 day period in blue is obtained from Sykulski et al posted Matlab code

#plot time series, coloring particular sections 
plot(drift$lon[startTs:endTs], drift$lat[startTs:endTs], type = "l", 
     xlab = "Longitude", ylab = "Latitude", main = "Replication of Figure 1 (Right)")
points(drift$lon[gStartTs:gEndTs], drift$lat[gStartTs:gEndTs], type = "l", col = "green")
points(drift$lon[bStartTs:bEndTs], drift$lat[bStartTs:bEndTs], type = "l", col = "blue")

###Figure 2###

#pull out complex valued velocities
bZ <- drift$cv[bStartTs:bEndTs]
gZ <- drift$cv[gStartTs:gEndTs]

#get periodogram
blue <- getPerio(Z = bZ, delta = DELTA)
green <- getPerio(Z = gZ, delta = DELTA)

#Find inertial frequency (PROBLEM HERE)
##TO DO: I AM CALCULATING THE INERTIAL FREQUENCY WRONG IN SOME WAY
#T <- (23*60*60) + (56*60) + 4.1 #time to rotate around Earth once
#psi <- 2*pi/T #rate of Earth's rotation
#bLat <- (pi/2)*(drift$lat[bStartTs:bEndTs]/90)
#bf0 <- 2*pi*-2*psi*sin(bLat)*1000
#mean(bf0)

#plot figure
par(mfrow = c(1, 2))
plot(DELTA*blue$omega, blue$sZ, type = "l", col = "blue", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi), ylim = c(-20, 60))
#abline(v = f0, col = "red")
plot(DELTA*green$omega, green$sZ, type = "l", col = "green", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi), ylim = c(-20, 60))
mtext("Replication of Figure 2", outer = T, line = -3)
#abline(v = f0, col = "red")

###Figure 3
###TO DO: How did they come up with these variance (autocorrelation at 0)?

#sample matern
N <- 1000
DELTA <- 2
matSamp <- simMatern(B = 10, alpha  = 0.9, h = 0.1, N = N, delta = DELTA)
mat <- getPerio(Z = matSamp, delta = DELTA)

#sample FBM
N <- 1000
fbmSamp <- simFBM(B = 10, alpha  = 0.9, H = 0.1, N = N)
fbm <- getPerio(fbmSamp, delta = DELTA)

#plot figure
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
plot(fbmSamp[,"u"], col = "blue", xlim = c(0, 1000), ylim = c(-200, 200),
     type = "l", xlab = expression(paste(t, Delta)), ylab = "u (cm/s)")
points(matSamp[,"u"], col = "red", type = "l")
mtext("Realization of Figure 3 \n (Note: Figure is not identical, since it represents a realization of a random processes)", 
      outer = T, line = -3)
plot(fbmSamp, col = "blue", xlim = c(-200, 200), ylim = c(-200, 200),
     xlab = "u (cm/s)", ylab = "v (cm/s)", type = "l")
points(matSamp, col = "red", type = "l")
plot(DELTA*fbm$omega, fbm$sZ, col = "blue", ylim = c(-20, 80), xlim = c(-pi, pi),
     type = "l", xlab = expression(paste(omega, Delta)), ylab = "dB")
points(DELTA*mat$omega, mat$sZ, col = "red", type = "l")

###Figure 4
CF <- mean(-4*pi*(drift$f[2971:3570]))
N <- length(bZ)
#Fit and plot blue time series
bFit <- fitModel(bZ, CF, delta = 2, fracNeg = 0.4, fracPos = 0, quantSet = .5)
par(mfrow = c(1,2))
plot(DELTA*blue$omega, blue$sZ, type = "l", col = "blue", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi), ylim = c(-20, 60))
mtext("Replication of Figure 4", outer = T, line = -3)
sTau <- ouAc(bFit$A, bFit$w0, bFit$C, N) + maternAc(bFit$B, bFit$alpha, bFit$h, N)  
tau <- seq(0, N - 1)
sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
points(DELTA*blue$omega, 10*log10(sBar), col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*blue$omega[bFit$firstIndex:bFit$lastIndex], 10*log10(sBar)[bFit$firstIndex:bFit$lastIndex], 
       col = "green", type = "l", lwd = 3)
abline(v = CF, lwd = 3) #inertial frequency
abline(v = bFit$w0, lwd = 3, lty = 3) #estimated w0
abline(v = CF - bFit$w0, lwd = 3, lty = 3) #shifted inertial freq (f = f0 - w_eddy)

#Fit and plot green time series
#abline(v = f0, col = "red")
gFit <- fitModel(gZ, CF, delta = 2, fracNeg = 0.4, fracPos = 0,  quantSet = 0.5)
plot(DELTA*green$omega, green$sZ, type = "l", col = "green", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi), ylim = c(-20, 60))
sTau <- ouAc(gFit$A, gFit$w0, gFit$C, N) + maternAc(gFit$B, gFit$alpha, gFit$h, N)  
tau <- seq(0, N - 1)
sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
points(DELTA*green$omega, 10*log10(sBar), col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*blue$omega[bFit$firstIndex:bFit$lastIndex], 10*log10(sBar)[bFit$firstIndex:bFit$lastIndex], 
       col = "blue", type = "l", lwd = 3)
abline(v = CF, lwd = 3) #inertial frequency

###Figure 5
newinertial <- readMat("/users/hdirector/Documents/prelim/prelim/Code/newInertial.mat")
DELTA <- 1
num = newinertial$newinertial[,,1]$drifters
vel <- num[[6]] + 1i*num[[7]] #complex velocities, measured in cm/s
vel1 <- vel[,1]

#inertial frequency
lat <- num[[10]] #latitudes
psi <- mean(lat[,1])*pi/180 #average lat in radians
f0  <- -(8*pi/23.9345)*sin(psi) #still need to figure out where these
                          #constant come from (physics, units issue)

#left figure (first of 200 simulated drifters)
sim1Fit <- fitModel(vel1, CF = f0, DELTA, 1.5*f0/pi, 1.5*f0/pi, quantSet = .9)  #fit model
sim1Per <- getPerio(vel1, delta = DELTA, dB = TRUE, noZero = FALSE) #periodogram
N <- length(vel1)
par(mfrow = c(1, 2))
plot(sim1Per$omega, sim1Per$sZ, type = "l", col = "blue",
     ylim = c(-40, 60), xlim = c(-pi, pi), ylab = "dB",
     xlab = expression(paste(omega, Delta)))
mtext("Replication of Figure 5", outer = T, line = -3)

sTau <- ouAc(sim1Fit$A, sim1Fit$w0, sim1Fit$C, N) + maternAc(sim1Fit$B, sim1Fit$alpha, sim1Fit$h, N)  
tau <- seq(0, N - 1)
sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
points(DELTA*sim1Per$omega, 10*log10(sBar), col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*sim1Per$omega[sim1Fit$firstIndex:sim1Fit$lastIndex], 10*log10(sBar)[sim1Fit$firstIndex:sim1Fit$lastIndex], 
       col = "green", type = "l", lwd = 3)

#right figure
#storage vector for 10log10(sBar), which we will need to calculate an average
nSim <- 200
sBarMat <- matrix(nrow = N, ncol = nSim)
sZMat <- matrix(nrow = N, ncol = nSim)
sBarMat[,1] <- 10*log10(sBar)
sZMat[,1] <- sim1Per$sZ

#Loop through and fit other 199 simulations 
for (i in 2:nSim) {
  #fit each of the 200 simulations
  tempFit <- fitModel(vel[,i], CF = f0, DELTA, 1.5*f0/pi, 1.5*f0/pi, quantSet = .9)  #fit model
  tempPer <- getPerio(vel[,i], delta = DELTA, dB = TRUE, noZero = FALSE) #periodogram
  N <- length(vel[,i])
  #plot spectrum and model fit 
  sTau <- ouAc(tempFit$A, tempFit$w0, tempFit$C, N) + maternAc(tempFit$B, tempFit$alpha, tempFit$h, N)  
  tau <- seq(0, N - 1)
  sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
  sBarMat[,i] <- 10*log10(sBar)
  sZMat[,i] <- tempPer$sZ
  print(i)
}

#plot first periodogram make initial plot
plot(sim1Per$omega, sim1Per$sZ, type = "l", col = "lightgrey",
     ylim = c(-40, 60), xlim = c(-pi, pi), ylab = "dB",
     xlab = expression(paste(omega, Delta)))

#Add all the other periodograms
for (i in 2:nSim) {
  points(sim1Per$omega, sZMat[,i], type = "l", col = "lightgray")
}
#plot all fits on top 
for (i in 1:nSim) {
  points(DELTA*sim1Per$omega, sBarMat[,i], col = "darkgrey", type = "l", lwd = 1)
}

#Calculate and add ensemble mean periodogram and fit
meanEnsemSz <- apply(sZMat, 1, mean)
meanEnsemSBar <- apply(sBarMat, 1, mean)
points(DELTA*sim1Per$omega, meanEnsemSz, col = "blue", type = "l")
points(DELTA*sim1Per$omega, meanEnsemSBar, col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*sim1Per$omega[sim1Fit$firstIndex:sim1Fit$lastIndex],
       meanEnsemSBar[sim1Fit$firstIndex:sim1Fit$lastIndex], col = "green", type = "l", lwd = 3)

###Figure 6
#12 measurements per day
drifterulysses <- readMat("/users/hdirector/Documents/prelim/prelim/Code/drifterulysses.mat")
driftUlys <- list("num" = as.vector(drifterulysses$drifterulysses[,,1]$num),
              "lat" = as.vector(drifterulysses$drifterulysses[,,1]$lat),
              "lon" = as.vector(drifterulysses$drifterulysses[,,1]$lon),
              "cv" = as.vector(drifterulysses$drifterulysses[,,1]$cv),
              "f" = as.vector(drifterulysses$drifterulysses[,,1]$f))
DELTA <- 2

par(mfrow = c(1, 2))
plot(driftUlys$lon, driftUlys$lat, type = "l", col = "blue", xlim = c(-120, -80),
     ylim = c(-40, -20), xlab = "Longitude", ylab = "Latitude") 
mtext("Replication of Figure 6", outer = T, line = -3)
#days of interest (day 358 to 378, 12 readings per day)
currIndex <- (350*12 + 1):(370*12) #matching code not paper (days 350 to 370, vs 358 to 378)
CF <- mean(4*pi*(driftUlys$f[currIndex]))
NSamp <- length(currIndex)
points(driftUlys$lon[currIndex], driftUlys$lat[currIndex], type = "l", col = "red")

#get Periodogram
samp <- getPerio(driftUlys$cv[currIndex], DELTA, dB = TRUE, noZero = TRUE)
sampFit <- fitModel(driftUlys$cv[currIndex], CF, delta = 2, fracNeg = 0, fracPos = .4,
                    quantSet = .5, incZero = TRUE)
plot(DELTA*samp$omega, samp$sZ, type = "l", col = "blue", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi), ylim = c(-10, 50))
sTau <- ouAc(sampFit$A, sampFit$w0, sampFit$C, NSamp) + maternAc(sampFit$B, sampFit$alpha, sampFit$h, NSamp)  
tau <- seq(0, NSamp - 1)
sBar <- 2*fft(sTau*(1 - (tau/NSamp))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
points(DELTA*samp$omega, 10*log10(sBar), col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*samp$omega[sampFit$firstIndex:sampFit$lastIndex], 10*log10(sBar)[sampFit$firstIndex:sampFit$lastIndex], 
       col = "green", type = "l", lwd = 3)
abline(v = CF, lwd = 3) #inertial frequency

#Figures 7-10
#constants
CF <- max(4*pi*(driftUlys$f)) #positive because of southern hemisphere, why at max not mean
N <- length(driftUlys$num)
DELTA <- 2
nWin <- 999

#find region to fit model in
fracPos <- 1.75*max(4*pi*driftUlys$f)/pi #figure out exactly the logic of this cut-off

#initial periodogram and par fit (use for set up)
toStartPer <- getPerio(driftUlys$cv[1:nWin], DELTA, dB = TRUE, noZero = TRUE) 
toStartFit <- fitModel(driftUlys$cv[1:nWin], CF,  DELTA, fracNeg = 0, fracPos = fracPos,
                       quantSet = 0.8, incZero = TRUE, hess = TRUE)

#frequencies of interest
firstIndex <- toStartFit$firstIndex
lastIndex <- toStartFit$lastIndex
nOmega <-length(firstIndex:lastIndex)
sZMat <- matrix(ncol = nOmega, nrow = N - 499*2)

#rolling window CF (taking the max)
CFVec <- rep(NA, N - 499*2)
for (i in 500:(N - 499)) {
  CFVec[i - 499] <- mean(4*pi*(driftUlys$f[(i - 499):(i + 499)])) #positive because of southern hemisphere, why at max not mean
}

#get observed periodogram
for (i in 500:(N - 499)) {
  tempPer <- getPerio(driftUlys$cv[(i - 499):(i + 499)], DELTA, dB = TRUE, noZero = TRUE) 
  sZMat[i - 499, ] <- tempPer$sZ[firstIndex:lastIndex]
}

#storage vectors for parameters and CI Bounds
par6Val <- matrix(nrow = N - 499*2, ncol = 6)
colnames(par6Val) <- c("A", "B", "w0", "C", "h", "alpha")
par6LB <- par6UB <- par6Val
par6Val[1, ] <- c(toStartFit$A, toStartFit$B, toStartFit$w0,
                  toStartFit$C, toStartFit$h, toStartFit$alpha)
#Confidence interval for first time point
step <- qnorm(.975)*sqrt(diag(solve(toStartFit$hess))/nWin)
par6LB[1, ] <- par6Val[1, ] - step; par6UB[1, ] <- par6Val[1, ] + step

#Loop through all time points and calc par estimate and CI
for (i in 501:(N - 499)) {
  currCv <- driftUlys$cv[(i - 499):(i + 499)] 
  #fit parameters, initialize estimates at value used in last run
  tempFit <- fitModel(currCv, CF, DELTA, fracNeg = 0, fracPos = fracPos,
                  quantSet = .5, needInits = FALSE, parInit = par6[i - 499 - 1, ],
                  hess = TRUE)
  par6Val[i - 499, ] <- c(tempFit$A, tempFit$B, tempFit$w0,
                    tempFit$C, tempFit$h, tempFit$alpha)
  #calculate confidence interval via Fisher's information
  step <- qnorm(.975)*sqrt(diag(solve(tempFit$hess))/nWin)
  par6LB[i - 499, ] <- par6Val[i - 499, ] - step; par6UB[i - 499, ] <- par6Val[i - 499, ] + step
  print(i)
}

#plot observed periodogram
par(mfrow = c(2, 1))
image.plot((500:(N - 499))/12, 2*toStartPer$omega[firstIndex:lastIndex],
           sZMat, useRaster = TRUE, zlim = c(20, 60),
           ylab = expression(paste("Frequency in radians (", omega, Delta %in% Omega, ")")),
           xlab = "Day")
points((500:(N - 499))/12, CFVec, col = "white", lwd = .01)
mtext("Replication of Figure 7", outer = T, line = -3)

#find and plot fitted periodogram
sZMatObsa <- matrix(ncol = nOmega, nrow = N - 499*2)
for (i in 1:(N - 499*2)) {
  sTau <- (ouAc(parVal[i, "A"], parVal[i, "w0"], parVal[i, "C"], nWin) +
                 maternAc(parVal[i, "B"], parVal[i, "alpha"], parVal[i, "h"], nWin))  
  tau <- seq(0, nWin - 1)
  sBar <- 2*fft(sTau*(1 - (tau/nWin))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
  sZMatObs[i, ] <- 10*log10(sBar)[toStartFit$firstIndex:toStartFit$lastIndex]
}

#Figure 8 parameter estimates over times
par(mfrow = c(1, 1))
plot((500:(N - 499))/12, par6[,"w0"], ylim = c(.43, .67), type = "l", 
     xlab = "Day", ylab = "Inertial Frequencies", col = "blue")
polygon(c(rev((500:(N - 499))/12), (500:(N - 499))/12), border = NA,
        c(rev(par6UB[, "w0"]), par6LB[, "w0"]), col = "lavender")
mtext("Replication of Figure 8", outer = T, line = -3)
points((500:(N - 499))/12, CFVec, col = "red", type ="l")
plot((500:(N - 499))/12, par6[, "A"], ylim = c(0, 27), type = "l",
     xlab = "Day", ylab  = "Amplitudes", col = "blue")
polygon(c(rev((500:(N - 499))/12), (500:(N - 499))/12),  border = NA,
        c(rev(par6UB[, "A"]), par6LB[, "A"]), col = "lavender")
points((500:(N - 499))/12, par6[, "B"], col = "red", type = "l")
polygon(c(rev((500:(N - 499))/12)[1:500], (500:(N - 499))/12)[1:500],  border = NA,
        c(rev(par6LB[, "B"])[1:500], par6UB[, "B"])[1:500], col = "pink")
plot((500:(N - 499))/12, par6[, "C"], ylim = c(0, 0.2), type = "l",
     xlab = "Day", ylab = "Dampening", col = "blue")
points((500:(N - 499))/12, par6[, "h"], col = "red", type = "l")
plot((500:(N - 499))/12, par6[, "alpha"], ylim = c(0.58, 1.45), type = "l",
     xlab = "Day", ylab = "Slope Parameter", col = "blue")



##Figure 9
#storage vectors for 5 parameter model
par5 <- matrix(nrow = N - 499*2, ncol = 5)
colnames(par5) <- c("A", "B", "w0", "C", "h", "alpha")