###Hannah Director###
#Research paper prelim 
rm(list = ls())
source("/users/hdirector/Documents/prelim/prelim/Code/functions.R")

###Figure 1 (left)###
#Image from global drifter project (paper author's did not make this artistic figure 
#or provide the data for it

############################################
#Figures 1-2 & 4, Nothern Hemisphere drifter
############################################

###Figure 1 (right)###
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

#Indices of sections of interest (b denotes the "blue" highlighted blue section and g the highlighted "green" section of the time series)
startTs <- 1901; endTs <- 4900 #index of period of interest is obtained from Sykulski et al. posted Matlab code
bStartTs <- 2971; bEndTs <-3570 #index of 50 day period in blue is obtained from Sykulski et al. posted Matlab code
gStartTs <- 4201; gEndTs <- 4800 #index of 50 day period in blue is obtained from Sykulski et al. posted Matlab code

#plot time series, coloring particular sections 
#pdf('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig1.pdf',
#    height = 3, width = 4)
#par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(5,5,4,1.1)) #for pretty exporting
par(mfrow = c(1, 1))
plot(drift$lon[startTs:endTs], drift$lat[startTs:endTs], type = "l", 
     xlab = "Longitude", ylab = "Latitude", main = "Replication of Figure 1 (Right)")
points(drift$lon[gStartTs:gEndTs], drift$lat[gStartTs:gEndTs], type = "l", col = "green")
points(drift$lon[bStartTs:bEndTs], drift$lat[bStartTs:bEndTs], type = "l", col = "blue")
legend("bottomright", legend = c("exhibits eddies", "no major eddies"), fill = c("blue", "green"), cex = .6)
#dev.off()

###Figure 2###
#Constants 
DELTA <- 2

#pull out complex valued velocities
bZ <- drift$cv[bStartTs:bEndTs]
gZ <- drift$cv[gStartTs:gEndTs]

#get periodogram
blue <- getPerio(Z = bZ, delta = DELTA)
green <- getPerio(Z = gZ, delta = DELTA)

#Find inertial frequency
gCf <- mean(-4*pi*drift$f[gStartTs:gEndTs])
bCf <- mean(-4*pi*drift$f[bStartTs:gStartTs])

#time series lengths
bN <- length(bZ)
gN <- length(gZ)

#plot figure
#pdf('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig2.pdf',
#    height = 3, width = 8.5)
#par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(4,5,4,1.1)) #for pretty exporting
par(mfrow = c(1, 2))
plot(DELTA*blue$omega, blue$sZ, type = "l", col = "blue", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi), ylim = c(-20, 60))
abline(v = bCf, col = "red")
legend("topright", legend = c("data", "inertial freq."),
       lty = c("solid", "solid"), col = c("blue", "red"),
       cex = .7, bg = "white")
plot(DELTA*green$omega, green$sZ, type = "l", col = "green", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi), ylim = c(-20, 60))
mtext("Replication of Figure 2", outer = T, line = -3)
abline(v = gCf, col = "red")
legend("topright", legend = c("data", "inertial freq."),
       lty = c("solid", "solid"), col = c("green", "red"),
       cex = .7, bg = "white")
#dev.off()

###Figure 4###
#get mean coriolis frequency for blue and green time series periods
bCF <- mean(-4*pi*(drift$f[bStartTs:bEndTs]))
gCF <- mean(-4*pi*(drift$f[gStartTs:gEndTs]))

#Fit and plot blue time series
bFit <- fitModel(bZ, bCF, delta = DELTA, fracNeg = 0.4, fracPos = 0, quantSet = .5)
#pdf('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig4.pdf',
#    height = 4, width = 8.5)
#par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(4,5,4,1.1)) #for pretty exporting
par(mfrow = c(1,2))
plot(DELTA*blue$omega, blue$sZ, type = "l", col = "blue", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi + 1), ylim = c(-20, 60))
mtext("Replication of Figure 4", outer = T, line = -3)
sTau <- ouAc(bFit$A, bFit$w0, bFit$C, bN, delta = DELTA) + maternAc(bFit$B, bFit$alpha, bFit$h, bN, delta = DELTA)  
tau <- seq(0, bN - 1)
sBar <- 2*fft(sTau*(1 - (tau/bN))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
points(DELTA*blue$omega, 10*log10(sBar), col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*blue$omega[bFit$firstIndex:bFit$lastIndex], 10*log10(sBar)[bFit$firstIndex:bFit$lastIndex], 
       col = "green", type = "l", lwd = 3)
abline(v = bCF, lwd = 3) #inertial frequency
abline(v = bFit$w0, lwd = 3, lty = 3, col = "gray") #estimated w0
abline(v = bCF - bFit$w0, lwd = 3, lty = 3) #shifted inertial freq (f = f0 - w_eddy)
legend("topright", legend = c("data", "fit", "extended fit", "eddy peak", "Coriolis freq.", "shifted inertial freq."),
       lty = c("solid", "solid", "dotted", "dotted","solid", "dotted"),
       col = c("blue", "green", "red", "gray", "black", "black"),
       cex = .7, bg = "white")

#Fit and plot green time series
gFit <- fitModel(gZ, gCF, delta = DELTA, fracNeg = 0.4, fracPos = 0,  quantSet = 0.5)
plot(DELTA*green$omega, green$sZ, type = "l", col = "green", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi + 1), ylim = c(-20, 60))
abline(v = gCF, col = "red")
sTau <- ouAc(gFit$A, gFit$w0, gFit$C, gN, delta = DELTA) + maternAc(gFit$B, gFit$alpha, gFit$h, gN, delta = DELTA)  
tau <- seq(0, gN - 1)
sBar <- 2*fft(sTau*(1 - (tau/gN))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
points(DELTA*green$omega, 10*log10(sBar), col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*blue$omega[bFit$firstIndex:bFit$lastIndex], 10*log10(sBar)[bFit$firstIndex:bFit$lastIndex], 
       col = "blue", type = "l", lwd = 3)
abline(v = gCF, lwd = 3) #inertial frequency
legend("topright", legend = c("data", "fit", "extended fit", "Coriolis freq."),
       lty = c("solid", "solid", "dotted","solid"),
       col = c("green", "blue", "red", "black"),
       cex = .7, bg = "white")
#dev.off()

##################################
#Figure 3: Simulate Matern and FBM
##################################
###Figure 3###
DELTA <- 2

#sample matern process
set.seed(103)
N <- 1000
DELTA <- 2
matSamp <- simMatern(B = 10, alpha  = 0.9, h = 0.1, N = N, delta = DELTA)

#sample fractional brownian motion process (fbm)
N <- 1000
fbmSamp <- simFBM(B = 10, alpha  = 0.9, H = 0.1, N = N)

#plot figure
pdf('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig3.pdf',
    height = 6, width = 6)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(2,2,3,1.1), mai = rep(.25, 4), oma = rep(1, 4)) #for pretty exporting
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
#u coordinate vs. time
plot(fbmSamp[,"u"], col = "blue", xlim = c(0, 1000), ylim = c(-200, 200),
     type = "l", xlab = expression(paste(t, Delta)), ylab = "u (cm/s)")
legend("topright", legend = c("FBM", "Matern"), fill = c("blue", "red"))
points(matSamp[,"u"], col = "red", type = "l")
mtext("Replication of Figure 3 (with a different realization)", 
      outer = T, line = -1)
#u coordinates vs. v coordinates
plot(fbmSamp, col = "blue", xlim = c(-200, 200), ylim = c(-200, 200),
     xlab = "u (cm/s)", ylab = "v (cm/s)", type = "l")
legend("topright", legend = c("FBM", "Matern"), fill = c("blue", "red"))
points(matSamp, col = "red", type = "l")
#periodogram
matSamp <- matSamp[,1] + matSamp[,2]*1i
mat <- getPerio(Z = matSamp, delta = DELTA)
fbmSamp <- fbmSamp[,1] + fbmSamp[,2]*1i
fbm <- getPerio(fbmSamp, delta = DELTA)
plot(DELTA*fbm$omega, fbm$sZ, col = "blue", ylim = c(-20, 80), xlim = c(-pi, pi),
     type = "l", xlab = expression(paste(omega, Delta)), ylab = "dB")
legend("topright", legend = c("FBM", "Matern"), fill = c("blue", "red"), cex = 1)
points(DELTA*mat$omega, mat$sZ, col = "red", type = "l")
dev.off()

#####################################
#Figure 5: SPEM Computer Experiment
#####################################
###Figure 5###
#read in data from SPEM computer experiment and add related info
newinertial <- readMat("/users/hdirector/Documents/prelim/prelim/Code/newInertial.mat")
num = newinertial$newinertial[,,1]$drifters
DELTA <- 1 #constant sampling increment
vel <- num[[6]] + 1i*num[[7]] #complex velocities, measured in cm/s
vel1 <- vel[,1] #first simulation, pulled out for left figure

#inertial frequency
lat <- num[[10]] 
psi <- mean(lat[,1])*pi/180 
f0  <- -(8*pi/23.9345)*sin(psi) 

#Consider only first of 200 simulated drifters for left figure
#Fit spectral model and plot periodogram and fitted model
sim1Fit <- fitModel(vel1, CF = f0, DELTA, 1.5*f0/pi, 1.5*f0/pi, quantSet = .8)  #fit model
sim1Per <- getPerio(vel1, delta = DELTA, dB = TRUE, noZero = FALSE) #periodogram
N <- length(vel1)
#png('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig5.png',
#    width = 8.5, height = 4, units = 'in', res = 300)
#par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(4,5,4,1.1)) #for pretty exporting
par(mfrow = c(1, 2))
plot(sim1Per$omega, sim1Per$sZ, type = "l", col = "blue",
     ylim = c(-40, 60), xlim = c(-pi, pi + 1), ylab = "dB",
     xlab = expression(paste(omega, Delta)))
mtext("Replication of Figure 5", outer = T, line = -3)
sTau <- ouAc(sim1Fit$A, sim1Fit$w0, sim1Fit$C, N, delta = DELTA) + maternAc(sim1Fit$B, sim1Fit$alpha, sim1Fit$h, N, delta = DELTA)  
tau <- seq(0, N - 1)
sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
points(DELTA*sim1Per$omega, 10*log10(sBar), col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*sim1Per$omega[sim1Fit$firstIndex:sim1Fit$lastIndex], 10*log10(sBar)[sim1Fit$firstIndex:sim1Fit$lastIndex], 
       col = "green", type = "l", lwd = 3)
legend("topright", legend = c("data", "fit", "extended fit"),
       lty = c("solid", "solid", "dotted"),
       col = c("blue", "green", "red"),
       cex = .7, bg = "white")

#right figure (all 200 simulations)
#storage vector for 10log10(sBar), from which we will need to calculate an average
nSim <- 200
sBarMat <- matrix(nrow = N, ncol = nSim)
sZMat <- matrix(nrow = N, ncol = nSim)

#store first value since we have already calculated it
sBarMat[,1] <- 10*log10(sBar)
sZMat[,1] <- sim1Per$sZ

#load results rather than re-run
load('/users/hdirector/Documents/prelim/prelim/Code/sZMat.rda')
load('/users/hdirector/Documents/prelim/prelim/Code/sBarMat.rda')

#Loop through other 199 simulations results and fit them
#for (i in 2:nSim) {
#  tempFit <- fitModel(vel[,i], CF = f0, DELTA, 1.5*f0/pi, 1.5*f0/pi, quantSet = .9)  
#  tempPer <- getPerio(vel[,i], delta = DELTA, dB = TRUE, noZero = FALSE) #periodogram
#  N <- length(vel[,i])
#  sTau <- ouAc(tempFit$A, tempFit$w0, tempFit$C, N, delta = DELTA) + maternAc(tempFit$B, tempFit$alpha, tempFit$h, N, delta = DELTA)  
#  tau <- seq(0, N - 1)
#  sBar <- 2*fft(sTau*(1 - (tau/N))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
#  sBarMat[,i] <- 10*log10(sBar)
#  sZMat[,i] <- tempPer$sZ
#  print(i)
#}

#store or load results
#save(sZMat, file = '/users/hdirector/Documents/prelim/prelim/Code/sZMat.rda')
#save(sBarMat, file = '/users/hdirector/Documents/prelim/prelim/Code/sBarMat.rda')


#plot first periodogram to make initial plot
plot(sim1Per$omega, sim1Per$sZ, type = "l", col = "lightgrey",
     ylim = c(-40, 60), xlim = c(-pi, pi + 1), ylab = "dB",
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
legend("topright", legend = c("periodograms", "mean periodogram", "fits", "mean fit", "mean fit extended"),
       lty = c("solid", "solid", "solid", "solid","dotted"),
       col = c("lightgrey", "blue", "darkgrey", "green", "red"),
       cex = .7, bg = "white")
#dev.off()

############################################
#Figures 6-10; Southern Hemisphere drifter
############################################
#Read in Southern Hemisphere drifter and related info data (has 12 measurements per day)
drifterulysses <- readMat("/users/hdirector/Documents/prelim/prelim/Code/drifterulysses.mat")
driftUlys <- list("num" = as.vector(drifterulysses$drifterulysses[,,1]$num),
                  "lat" = as.vector(drifterulysses$drifterulysses[,,1]$lat),
                  "lon" = as.vector(drifterulysses$drifterulysses[,,1]$lon),
                  "cv" = as.vector(drifterulysses$drifterulysses[,,1]$cv),
                  "f" = as.vector(drifterulysses$drifterulysses[,,1]$f))
#constants
DELTA <- 2 
MPD <- 12 #measurements per day
N <- length(driftUlys$num)
nWin <- 1000 #rolling window length
xVal <- (500:(N - 500))/MPD #Convert x-axis to days not number of observations

###Figure 6###
#pdf('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig6.pdf',
#    height = 3, width = 8.5)
#par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(4,5,4,1.1)) #for pretty exporting
par(mfrow = c(1, 2))
#plot langrangian time series
plot(driftUlys$lon, driftUlys$lat, type = "l", col = "blue", xlim = c(-120, -80),
     ylim = c(-40, -20), xlab = "Longitude", ylab = "Latitude") 
mtext("Replication of Figure 6", outer = T, line = -3)
legend("top", legend = c("1642-day trajectory", "Days 350-370"),
       lty = c("solid", "solid"), col = c("blue", "red"),
       cex = .7, bg = "white", horiz = TRUE)

#Specify days of interest, Note: I use days 350 to 370 not days 358 to 378. This matches the 
#posted code and figure in the paper, not the figure caption 
currIndex <- (350*12 + 1):(370*12) 
CF <- mean(4*pi*(driftUlys$f[currIndex]))
NSamp <- length(currIndex)
points(driftUlys$lon[currIndex], driftUlys$lat[currIndex], type = "l", col = "red")

#Plot periodogram and model period for highlighted time period of interest
samp <- getPerio(driftUlys$cv[currIndex], DELTA, dB = TRUE, noZero = TRUE)
sampFit <- fitModel(driftUlys$cv[currIndex], CF, delta = 2, fracNeg = 0, fracPos = .4,
                    quantSet = .5, incZero = FALSE, getHess = TRUE)
plot(DELTA*samp$omega, samp$sZ, type = "l", col = "blue", xlab = expression(paste(omega, Delta)),
     ylab = "dB", xlim = c(-pi, pi), ylim = c(-10, 50))
sTau <- ouAc(sampFit$A, sampFit$w0, sampFit$C, NSamp, delta = DELTA) + maternAc(sampFit$B, sampFit$alpha, sampFit$h, NSamp, delta = DELTA)  
tau <- seq(0, NSamp - 1)
sBar <- 2*fft(sTau*(1 - (tau/NSamp))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
points(DELTA*samp$omega, 10*log10(sBar), col = "red", lty = 3, type = "l", lwd = 3)
points(DELTA*samp$omega[sampFit$firstIndex:sampFit$lastIndex], 10*log10(sBar)[sampFit$firstIndex:sampFit$lastIndex], 
       col = "green", type = "l", lwd = 3)
abline(v = CF, lwd = 3) #inertial frequency
legend("topright", legend = c("data", "fit", "extended fit"),
       lty = c("solid", "solid", "dotted"), col = c("blue", "green", "red"),
       cex = .7, bg = "white", horiz = FALSE)
#dev.off()

###Fit 6-parameter model, used throughout sectiont###
#find region to fit model in (semi-parametric, choice based on understanding of the physics)
fracPos <- 1.75*max(4*pi*driftUlys$f)/pi 

#Periodogram and paremeter fit for first rolling window (used for set up)
toStartPer <- getPerio(driftUlys$cv[1:nWin], DELTA, dB = TRUE, noZero = TRUE) 
toStartFit <- fitModel(driftUlys$cv[1:nWin], CF,  DELTA, fracNeg = 0, fracPos = fracPos,
                       quantSet = 0.8, incZero = FALSE, getHess = TRUE)
step <- qnorm(.975)*sqrt(diag(solve(toStartFit$hess)))
firstIndex <- toStartFit$firstIndex
lastIndex <- toStartFit$lastIndex
nOmega <-length(firstIndex:lastIndex)
sZMat <- matrix(ncol = nOmega, nrow = N - 999)

#Calculate the mean coriolis frequency for each moving window
CFVec <- rep(NA, N - 999)
for (i in 500:(N - 500)) {
  CFVec[i - 499] <- mean(4*pi*(driftUlys$f[(i - 499):(i + 500)])) #positive because of southern hemisphere
}

#get observed periodogram for each moving window (only include frequencies which are used in fitting)
for (i in 500:(N - 500)) {
  currCv <- driftUlys$cv[(i - 499):(i + 500)] 
  tempPer <- getPerio(currCv, 2, dB = TRUE, noZero = TRUE)
  sZMat[i - 499, ] <- tempPer$sZ[firstIndex:lastIndex]
}
sZMat[sZMat < 20] <- 20

#Storage vectors for parameters, negative log likelihood, confidence interval bands, and hessians
par6Val <- matrix(nrow = N - 999, ncol = 6)
colnames(par6Val) <- c("A", "B", "w0", "C", "h", "alpha")
par6LB <- par6UB <- par6Val
llVal6 <- rep(NA, N - 999)
llVal6 <- toStartFit$llVal
hessArray <- array(dim = c(N - 999, 6, 6), data = NA)
hessArray[1,,] <- toStartFit$hess
par6LB[1, ] <- par6Val[1, ] - step; par6UB[1, ] <- par6Val[1, ] + step

#store data for first window 
par6Val[1, ] <- c(toStartFit$A, toStartFit$B, toStartFit$w0,
                  toStartFit$C, toStartFit$h, toStartFit$alpha)
par6LB[1, ] <- par6Val[1, ] - step
par6UB[1, ] <- par6Val[1, ] + step

#load results, rather than re-running
load('/users/hdirector/Documents/prelim/prelim/Code/par6Val.rda')
load('/users/hdirector/Documents/prelim/prelim/Code/par6LB.rda')
load('/users/hdirector/Documents/prelim/prelim/Code/par6UB.rda')
load('/users/hdirector/Documents/prelim/prelim/Code/llVal6.rda')
load('/users/hdirector/Documents/prelim/prelim/Code/hess.rda')

#Loop through all windows and calculate parameter estimates, confidence intervals, and hessian
#for (i in 501:(N - 500)) {
#  currCv <- driftUlys$cv[(i - 499):(i + 500)] 
#  CFCurr <- CFVec[i - 499]
#fit parameters, initialize estimates at value used in last run
#  tempFit <- fitModel(currCv, CFCurr, DELTA, fracNeg = 0, fracPos = fracPos,
#                      quantSet = .5, needInits = FALSE, parInit = par6Val[i - 499 - 1, ],
#                      getHess = TRUE)
#  par6Val[i - 499, ] <- c(tempFit$A, tempFit$B, tempFit$w0,
#                          tempFit$C, tempFit$h, tempFit$alpha)
#  llVal6[i - 499] <-  tempFit$llVal
#calculate confidence interval via Fisher's information
#  step <- qnorm(.975)*sqrt(diag(solve(tempFit$hess)))
#  par6LB[i - 499, ] <- par6Val[i - 499, ] - step; par6UB[i - 499, ] <- par6Val[i - 499, ] + step
#  hessArray[i - 499,,] <- tempFit$hess
#  print(i)
#}

#store or load results
#save(par6Val, file = '/users/hdirector/Documents/prelim/prelim/Code/par6Val.rda')
#save(par6LB, file =  '/users/hdirector/Documents/prelim/prelim/Code/par6LB.rda')
#save(par6UB, file = '/users/hdirector/Documents/prelim/prelim/Code/par6UB.rda')
#save(llVal6, file = '/users/hdirector/Documents/prelim/prelim/Code/llVal6.rda')
#save(hessArray, file = '/users/hdirector/Documents/prelim/prelim/Code/hess.rda')

#plot observed periodogram
#png('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig7.png',
#  width = 8.5, height = 5, units = 'in', res = 300)
#par(mgp=c(.95, 0.45,0), mar=c(1, 3,2, 1), mai = c(.4, .85, .32, .25), oma = rep(.5, 4))
par(mfrow = c(2, 1))

library("fields")
image.plot(xVal, 2*toStartPer$omega[firstIndex:lastIndex],
           sZMat, useRaster = TRUE, zlim = c(20, 60),
           ylab = expression(paste("Freq. in radians (", omega, Delta %in% Omega, ")")),
           xlab = "Day", cex.lab = 1, cex.axis = .1)
points(xVal, CFVec, col = "black", lwd = 0.005)
mtext("Replication of Figure 7", outer = T, line = -1.5)

#find fitted periodogram
sZMatObs <- matrix(ncol = nOmega, nrow = N - 999)
for (i in 1:(N - 999)) {
  sTau <- (ouAc(par6Val[i, "A"], par6Val[i, "w0"], par6Val[i, "C"], nWin, delta = DELTA) +
             maternAc(par6Val[i, "B"], par6Val[i, "alpha"],  par6Val[i,"h"], nWin, delta = DELTA))  
  tau <- seq(0, nWin - 1)
  sBar <- 2*fft(sTau*(1 - (tau/nWin))) - sTau[1]; sBar = abs(Re(fftshift(sBar))) 
  sZMatObs[i, ] <- 10*log10(sBar)[toStartFit$firstIndex:toStartFit$lastIndex]
}

#plot fitted periodogram
image.plot(xVal, 2*toStartPer$omega[firstIndex:lastIndex],
           sZMatObs, useRaster = TRUE, zlim = c(20, 60),
           ylab = expression(paste("Freq. in radians (", omega, Delta %in% Omega, ")")),
           xlab = "Day", cex.lab = 1, cex.axis = 1)
points(xVal, CFVec, col = "black", lwd = 0.005)
#dev.off()

###Figure 8###
#plot Figures
#png('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig8.png',
#    width = 8.5, height = 8, units = 'in', res = 300)
#par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.4,5,4,1.1)) #for pretty exporting
par(mfrow = c(4, 1))

#plot w0 over time with confidence interval and coriolis frequency over time
plot(xVal, par6Val[,"w0"], ylim = c(.43, .67), type = "l", cex.axis = 2,
     xlab = "Day", ylab = "Inertial Freq.", col = "blue", cex.lab = 2) #initialize plot
polygon(c(rev(xVal), xVal), c(rev(par6LB[, "w0"]), par6UB[, "w0"]), col = 'lightblue', border = NA)
points(xVal, par6Val[,"w0"],  type = "l", col = "blue") #plot again, since line gets covered by polygon
points(xVal, CFVec, col = "red", type ="l")
mtext("Replication of Figure 8", outer = T, line = -2,
      cex = 1.5)
legend("top", legend = c(expression(hat(omega)[0]), "95% CI" , "Coriolis freq. "),
       fill = c("blue", "lightblue", "red"), horiz = TRUE, cex = 1.4)

#plot inertial frequencies over time with confidence interval 
plot(xVal, par6Val[, "A"], ylim = c(0, 32), type = "l", cex.axis = 2,
     xlab = "Day", ylab  = "Amplitudes", col = "blue", cex.lab = 2)
polygon(c(rev(xVal), xVal), c(rev(par6LB[, "A"]), par6UB[, "A"]), col = 'lightblue', border = NA)
points(xVal, par6Val[, "A"], type = "l", col = "blue")
polygon(c(rev(xVal), xVal), c(rev(par6LB[, "B"]), par6UB[, "B"]), col = 'pink', border = NA)
points(xVal, par6Val[, "B"], col = "red", type = "l")
legend("top", legend = c(expression(hat(A)),"95% CI", expression(hat(B)), "95% CI"), cex = 1.4,
       fill = c("blue", "lightblue", "red", "lightpink"), horiz = TRUE)

#plot dampening over time with confidence interval
plot(xVal, par6Val[, "h"], ylim = c(0, 0.16), type = "l", cex.axis = 2,
     xlab = "Day", ylab = "Dampening", col = "red", cex.lab = 2)
polygon(c(rev(xVal), xVal), c(rev(par6LB[, "h"]), par6UB[, "h"]), col = 'pink', border = NA)
polygon(c(rev(xVal), xVal), c(rev(par6LB[, "C"]), par6UB[, "C"]), col = 'lightblue', border = NA)
points(xVal, par6Val[, "C"], col = "blue", type = "l")
points(xVal, par6Val[, "h"], type = "l", col = "red")
legend("top", legend = c(expression(hat(c)),"95% CI", expression(hat(h)), "95% CI"), cex = 1.4,
       fill = c("blue", "lightblue", "red", "lightpink"), horiz = TRUE)

#plot slope over time with confidence interval
plot(xVal, par6Val[, "alpha"], ylim = c(0.6, 1.5), type = "l", cex.axis = 2, 
     xlab = "Day", ylab = "Slope", col = "blue", cex.lab = 2)
polygon(c(rev(xVal), xVal), c(rev(par6LB[, "alpha"]), par6UB[, "alpha"]), col = 'lightblue', border = NA)
points(xVal, par6Val[, "alpha"],  type = "l",col = "blue")
legend("top", legend = c(expression(hat(alpha)), "95% CI"), cex = 1.4,
       fill = c("blue", "lightblue"), horiz = TRUE)
#dev.off()

###Figure 9 ####
#load results, rather than re-running
#load('/users/hdirector/Documents/prelim/prelim/Code/llVal5.rda')

#Fraction positive 
fracPos <- 1.75*max(4*pi*driftUlys$f)/pi 

#storage vectors for negative log likelihood and parameters
N <- length(driftUlys$cv)
llVal5 <- rep(NA, N - 999)
par5Val <- matrix(nrow = N - 499, ncol = 5, data = NA)

#fit and store results for first time window (for set up)
toStartFit5 <- fitModel(driftUlys$cv[1:1000], CFVec[1], DELTA, fracNeg = 0, fracPos = fracPos,
                        quantSet = .8, simpModel = TRUE)
llVal5[1] <- toStartFit5$llVal
par5Val[1, ] <- c(toStartFit5$A, toStartFit5$B, toStartFit5$C, toStartFit5$h, toStartFit5$alpha)

#load results, rather than re-run
load('/users/hdirector/Documents/prelim/prelim/Code/llVal5.rda')

#Loop through all time points, fit simplified model, and store only likelihood
#for (i in 501:(N - 999)) {
#  currCv <- driftUlys$cv[(i - 499):(i + 500)] 
#  CFCurr <- CFVec[i - 499]
#fit parameters, initialize estimates at value used in last run
#  tempFit <- fitModel(currCv, CFCurr, DELTA, fracNeg = 0, fracPos = fracPos,
#                      quantSet = .8, simpModel = TRUE, needInits = FALSE, parInit = par5Val[i - 499 -1, ])
#store negative log likelihood value and parameter 
#  llVal5[i - 499] <-  tempFit$llVal
#  par5Val[i - 499, ] <- c(toStartFit$A, toStartFit$B, toStartFit$C, toStartFit$h, toStartFit$alpha)
#  print(i)
#}

#save results
#save(llVal5, file = '/users/hdirector/Documents/prelim/prelim/Code/llVal5.rda')

#plot likelihood ratio test statistic
#pdf('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig9.pdf',
#    height = 2.5, width = 7)
#par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.4,5,4,1.1)) #for pretty exporting
par(mfrow = c(1, 1))
#Note that we have stored the negative of the log likelihoods
LRT <- 2*(-llVal6 - -llVal5)
plot(xVal, LRT, type = "l", col = "blue", xlab = "Day",
     ylab = "Likelihood Ratio \n Test Statistic", main = "Replication of Figure 9")
abline(h = qchisq(.95, 1), col = "red", lty = 2)
#dev.off()

###Figure 10###
#Get Fisher info from Hessian
fishArray <- array(dim = dim(hessArray), data = NA)
for (i in 1:dim(hessArray)[1]) {
  fishArray[i,,] <- cov2cor(solve(hessArray[i,,]))
}

#average over all time periods
medCorr <- apply(fishArray, c(2, 3), function(x){mean(x, na.rm = T)})
medCorr <- round(medCorr, 2)
colnames(medCorr) <- c("A", "B", "w0", "c", "h", "alpha")
rownames(medCorr) <- c("A", "B", "w0", "c", "h", "alpha")

#re-arrange to match paper
newCorr <- matrix(nrow = 6, ncol = 6, data = NA)
colnames(newCorr) <- c("A", "w0", "c", "B", "alpha", "h")
rownames(newCorr) <- c("A", "w0", "c", "B", "alpha", "h")
newCorr[1, 1] <- medCorr[2, 2] 
newCorr[1, 2] <- newCorr[2, 1] <- medCorr[1, 3]
newCorr[1, 3] <- newCorr[3, 1] <- medCorr[1, 4]
newCorr[1, 4] <- newCorr[4, 1] <- medCorr[1, 2]
newCorr[1, 5] <- newCorr[5, 1] <- medCorr[1, 6]
newCorr[1, 6] <- newCorr[6, 1] <- medCorr[1, 5]
newCorr[2, 2] <- medCorr[3, 3]
newCorr[2, 3] <-  newCorr[3, 2] <- medCorr[3, 4]
newCorr[2, 4] <- newCorr[4, 2] <- medCorr[3, 2]
newCorr[2, 5] <- newCorr[5, 2] <- medCorr[3, 6]
newCorr[2, 6] <- newCorr[6, 2] <- medCorr[3, 5]
newCorr[3, 3] <- medCorr[4, 4]
newCorr[3, 4] <- newCorr[4, 3] <- medCorr[2, 4]
newCorr[3, 5] <- newCorr[5, 3] <- medCorr[4, 6]
newCorr[3, 6] <- newCorr[6, 3] <- medCorr[4, 5]
newCorr[4, 4] <- medCorr[2, 2]
newCorr[4, 5] <- newCorr[5, 4] <- medCorr[2, 6]
newCorr[4, 6] <- newCorr[6, 4] <- medCorr[2, 5]
newCorr[5, 5] <- medCorr[6, 6]
newCorr[5, 6] <- newCorr[6, 5] <- medCorr[6, 5]
newCorr[6, 6] <- medCorr[5, 5]


#plot correlation
library(corrplot)
pdf('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig10.pdf')
corrplot(newCorr, method = "color", p.mat = newCorr, insig = "p-value",
         sig.level = -1, type = "lower")
mtext("Replication of Figure 10", outer = T, line = -3, cex = 2)
dev.off()

#####################################
#Extension plots, change window size
#####################################
###Make LRT plot with window length of 500####
#find region to fit model in (semi-parametric, choice based on understanding of the physics)
fracPos <- 1.75*max(4*pi*driftUlys$f)/pi 
N <- length(driftUlys$cv)

#Calculate the mean coriolis frequency for each moving window
CFVec <- rep(NA, N - 1999)
for (i in 999:(N - 1000)) {
  CFVec[i - 999] <- mean(4*pi*(driftUlys$f[(i - 999):(i + 1000)])) #positive because of southern hemisphere
}

#Periodogram and paremeter fit for first rolling window (used for set up)
toStartPer <- getPerio(driftUlys$cv[1:2000], DELTA, dB = TRUE, noZero = TRUE) 
toStartFit <- fitModel(driftUlys$cv[1:2000], mean(CFVec[1]),  DELTA, fracNeg = 0, fracPos = fracPos,
                       quantSet = 0.8, incZero = FALSE, getHess = TRUE)
step <- qnorm(.975)*sqrt(diag(solve(toStartFit$hess)))
firstIndex <- toStartFit$firstIndex
lastIndex <- toStartFit$lastIndex

#Storage vectors for parameters, negative log likelihood, confidence interval bands, and hessians
par6Val_2000 <- matrix(nrow = N - 1999, ncol = 6)
colnames(par6Val_500) <- c("A", "B", "w0", "C", "h", "alpha")
llVal6_2000 <- rep(NA, N - 1999)
llVal6_2000[1] <- toStartFit$llVal

#store data for first window 
par6Val_2000[1, ] <- c(toStartFit$A, toStartFit$B, toStartFit$w0,
                  toStartFit$C, toStartFit$h, toStartFit$alpha)

#Loop through all windows and calculate parameter estimates, confidence intervals, and hessian
for (i in 1001:(N - 1000)) {
  currCv <- driftUlys$cv[(i - 999):(i + 1000)] 
  CFCurr <- CFVec[i - 999]
#fit parameters, initialize estimates at value used in last run
  tempFit <- fitModel(currCv, CFCurr, DELTA, fracNeg = 0, fracPos = fracPos,
                      quantSet = .8, needInits = FALSE, parInit = par6Val_2000[i - 999 - 1, ],
                      getHess = FALSE)
  par6Val_2000[i - 999, ] <- c(tempFit$A, tempFit$B, tempFit$w0,
                          tempFit$C, tempFit$h, tempFit$alpha)
  llVal6_2000[i - 999] <-  tempFit$llVal
  print(i)
}

#store  results
save(llVal6_2000, file = '/users/hdirector/Documents/prelim/prelim/Code/llVal6_2000.rda')

###Simplified model using longer window
#Periodogram and paremeter fit for first rolling window (used for set up)
toStartPer <- getPerio(driftUlys$cv[1:2000], DELTA, dB = TRUE, noZero = TRUE) 
toStartFit <- fitModel(driftUlys$cv[1:2000], mean(CFVec[1]),  DELTA, fracNeg = 0, fracPos = fracPos,
                       quantSet = 0.8, incZero = FALSE, getHess = TRUE, simpModel = TRUE)
firstIndex <- toStartFit$firstIndex
lastIndex <- toStartFit$lastIndex

#Storage vectors for parameters, negative log likelihood, confidence interval bands, and hessians
par5Val_2000 <- matrix(nrow = N - 1999, ncol = 5)
colnames(par5Val_2000) <- c("A", "B", "C", "h", "alpha")
llVal5_2000 <- rep(NA, N - 1999)
llVal5_2000[1] <- toStartFit$llVal

#store data for first window 
par5Val_2000[1, ] <- c(toStartFit$A, toStartFit$B, 
                       toStartFit$C, toStartFit$h, toStartFit$alpha)

#Loop through all windows and calculate parameter estimates, confidence intervals, and hessian
#for (i in 1001:(N - 1000)) {
#  currCv <- driftUlys$cv[(i - 999):(i + 1000)] 
#  CFCurr <- CFVec[i - 999]
  #fit parameters, initialize estimates at value used in last run
#  tempFit <- fitModel(currCv, CFCurr, DELTA, fracNeg = 0, fracPos = fracPos,
#                      quantSet = .8, needInits = FALSE, parInit = par5Val_2000[i - 999 - 1, ],
#                      getHess = FALSE, simpModel = TRUE)
#  par5Val_2000[i - 999, ] <- c(tempFit$A, tempFit$B,
#                               tempFit$C, tempFit$h, tempFit$alpha)
#  llVal5_2000[i - 999] <-  tempFit$llVal
#  print(i)
#}

#save(llVal5_2000, file = '/users/hdirector/Documents/prelim/prelim/Code/llVal5_2000.rda')

#both plots
#Note that we have stored the negative of the log likelihoods
load('/users/hdirector/Documents/prelim/prelim/Code/llVal6.rda')
load('/users/hdirector/Documents/prelim/prelim/Code/llVal5.rda')
load('/users/hdirector/Documents/prelim/prelim/Code/llVal5_2000.rda')
load('/users/hdirector/Documents/prelim/prelim/Code/llVal6_2000.rda')
pdf('/users/hdirector/Documents/prelim/prelim/Paper/ReplicatedFigures/fig9Extend.pdf',
    height = 3, width = 4)
par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(1.6,3.5,1,1)) #for pretty exporting
par(mfrow = c(2, 1))
LRT <- 2*(-llVal6 - -llVal5)
plot(xVal, LRT, type = "l", col = "blue", xlab = "Day", ylim = c(0, 40),
     ylab = "Likelihood Ratio \n Test Statistic", 
     main = "Replication of Figure 9, Window Length: 1000", cex.main = .6,
     cex.lab = .6, cex.axis = .6)
abline(h = qchisq(.95, 1), col = "red", lty = 2)
legend("topright", legend = c("LRT Test Stat.", "Significance Level"), fill = c("blue", "red"), cex = .6)


LRT2 <- 2*(-llVal6_2000 - -llVal5_2000)
LRT2 <- c(rep(NA, 500), LRT2, rep(NA, 500)) #don't plot in window that we can't observe
plot(xVal, LRT2, type = "l", col = "blue", xlab = "Day", ylim = c(0, 40),
     ylab = "Likelihood Ratio \n Test Statistic", 
     main = "Analogous Plot to Figure 9, Window Length: 2000", cex.main = .6, 
     cex.lab = .6, cex.axis = .6)
abline(h = qchisq(.95, 1), col = "red", lty = 2)
legend("topright", legend = c("LRT Test Stat.", "Significance Level"), fill = c("blue", "red"), cex = .6)

dev.off()
