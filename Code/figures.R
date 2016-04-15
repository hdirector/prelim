###Hannah Director###
#Research paper prelim 
rm(list = ls())
source("/users/hdirector/Dropbox/Prelim/ReplicationCode/functions.R")

#Constants
DELTA <- 2

#Convert data from matlab to R
library("R.matlab")
drifterbetty2 <- readMat("/users/hdirector/Dropbox/Prelim/ReplicationCode/drifterbetty2.mat")
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
#To do: figure out how to get complex valued velocities (e.g. units)
#Calculate velocities u(t) (horizontal) and v(t) (vertical)

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
matSamp <- simMatern(B = 10, alpha  = 0.9, h = 0.1, N = N)
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
optim(par = c(5, 2, 2, 4, 2, 2), ll, Z = bZ, delta = 2,  method = "L-BFGS-B", control=list(fnscale = -1),
      lower = c(0, -Inf, 0, 0, 0.5, 0), upper = c(Inf, Inf, Inf, Inf, Inf, Inf))


