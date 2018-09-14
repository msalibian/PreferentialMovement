####Load the necessary libraries
library(RandomFields)
library(geoR)
library(fields)
library(mvtnorm)
library(numDeriv)
library(MASS)
library(foreach)
library(doParallel)
library(matrixcalc)
library(prodlim)
library(TMB)
library(INLA)
#################################################################################
# genPrefTrackHybrid ############################################################
#################################################################################
#alternative track generation based on gradient of temperature field, not the temperature itself.
genPrefTrackHybridBehav <- function(n, param, paramGP, nrowcol, m, surface, 
                                    gridFull, model, rawDat, 
                                    alpha, start, burnIn, timing){
  # initialise beta (behavioural state)
  betaF <- NULL
  betaF <- c(betaF, alpha[1])
  behavState <- exp(alpha[1])/(1+exp(alpha[1]))
  # measurement error
  sigmasq.nugget <- paramGP[3]
  # create vector for true surface
  S.vec =  as.numeric(surface)
  # mux and muy are x and y partial derivatives of H respectively
  # combined these are the drift function in the SDE. Not sure why this
  # function exists since its only valid for t=0 UPDATE
  mux <- function(x, y, param, tempGrad, temp, behavState){
    Hxy <- behavState*alpha[2]*as.numeric(tempGrad[1])*as.numeric(temp)
    return(Hxy) 
  }
  muy <- function(x, y, param, tempGrad, temp, behavState){
    Hxy <- behavState*alpha[3]*as.numeric(tempGrad[2])*as.numeric(temp)
    return(Hxy) 
  }
  # mu is the drift function in the SDE for first time point
  mu <- function(x, y, param, tempGrad, temp, behavState){
    return(cbind(mux(x, y, param, tempGrad, temp, behavState),muy(x, y, param, tempGrad, temp, behavState)))
  }
  # muNew is drift function for all points after the first time point
  muNew <- function(path, param, tempGrad, alpha, behavState){
    # x drift direction
    Hxy1 <- alpha[2]*as.numeric(tempGrad[1])*as.numeric(path[nrow(path),4])
    # y drift direction
    Hxy2 <- alpha[3]*as.numeric(tempGrad[2])*as.numeric(path[nrow(path),4])
    # find time increment 
    timeDiff <- path[nrow(path),1] - path[(nrow(path)-1),1]
    # estimate velocity
    diff <- c((path[nrow(path),2]-path[(nrow(path)-1),2])/timeDiff, 
              (path[nrow(path),3]-path[(nrow(path)-1),3])/timeDiff)
    # caluclate full drift function
    resx <- ((1-behavState)*diff[1] + behavState*Hxy1)
    resy <- ((1-behavState)*diff[2] + behavState*Hxy2)
    return(cbind(resx, resy))
  }
  # initialise the path
  path <- NULL
  # initalise both time and location
  # column 1 is time
  # column 2 is longitude (x)
  # column 3 is latitude (y)
  # column 4 is S (ie/ temperature)
  # column 5 is x gradient
  # column 6 is y gradient
  if(start==0){
    obs <- c(0,runif(1,-150,150),runif(1,-150,150))
  }else{
    obs <- c(0,rnorm(1,0,start), rnorm(1,0, start))
  }
  # calculate gradient of field at current location
  tempGrad <- tempDrift(xy=as.numeric(obs[2:3]), paramGP=paramGP, surface=S.vec, frame=gridFull, 
                        nrowcol=nrowcol)
  # calculate temperature at current location
  obsT <- findTempGPV(obs[2:3], paramGP, S.vec, gridFull[1:(nrowcol^2),], nrowcol, neighbours=50)
  # create matrix of locations and temperatures
  gridFullL <- rbind(gridFull, obs[2:3])
  S.vec <- c(S.vec, obsT)
  # update path to include temperature
  obs <- c(obs,obsT)
  time <- 0
  path <- rbind(path,c(obs, NA, NA))
  # simulate first movement
  repeat{ 
    # simulate time difference
    timePart <- rexp(1,timing) + 0.003
    # simulate change in behaviour state
    alpha[1] <- rnorm(1, alpha[1], param[1]*sqrt(timePart))
    behavState <- exp(alpha[1])/(1+exp(alpha[1]))
    # calculate movement error
    errX <- sqrt(timePart)*param[2]*mvrnorm(n = 1, c(0,0), diag(2), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    # update total elapsed time
    time <- time + timePart
    # update path of trajectory using Brillinger's SDE linear approximation
    obs <- c(time, as.numeric(path[1,2:3]) + (as.numeric(path[1,4])*mu(path[1,2], path[1,3], param, tempGrad, path[1,4], behavState)*
                                                (timePart)) + errX, path[1,4]) 
    if(abs(obs[2])<150 && abs(obs[3])<150) {
      break
    }
  }
  # update behavioural state vector
  betaF <- c(betaF, alpha[1])
  obsT <- findTempGPV(as.numeric(obs[2:3]), paramGP, S.vec, gridFullL[1:(nrowcol^2),], nrowcol, neighbours=50)
  tempGrad <- tempDrift(xy=as.numeric(obs[2:3]), paramGP=paramGP, surface=S.vec, frame=gridFullL[1:(nrowcol^2),], 
                          nrowcol=nrowcol)
  # include nugget error on temperature
  obs[4] <- rnorm(1, mean=obsT, sd=sqrt(sigmasq.nugget))
  path <- rbind(path,c(obs, tempGrad))
  gridFullL <- rbind(gridFullL, obs[2:3])
  S.vec <- c(S.vec, obsT)
  # now simulate track   
  # t=1 for non-equal time intervals
  # time intervals are currently drawn from iid exponential(1)
  for(i in 2:n){ 
      # simulate time step
      timePart <- rexp(1,timing) + 0.003
      # simulate behavioural state
      alpha[1] <- rnorm(1, alpha[1], param[1]*sqrt(timePart))
      behavState <- exp(alpha[1])/(1+exp(alpha[1]))
      # movement diffusion 
      errX <- sqrt(timePart)*param[2]*mvrnorm(n = 1, c(0,0), diag(2), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
      # update total elapsed time
      time <- time + timePart
      # update path of trajectory using Brillinger's SDE linear approximation
      obs <- c(time, as.numeric(path[i,2:3]) + (muNew(path, param, tempGrad, alpha, behavState)*(timePart)) + errX, path[i,4]) 
      # break if outside 'safezone'
      if(abs(obs[2])>150|abs(obs[3])>150){
        break
      }else{
          obsT <- findTempGPV(as.numeric(obs[2:3]), paramGP, S.vec, gridFullL[1:(nrowcol^2),], nrowcol, neighbours=50)
          tempGrad <- tempDrift(xy=as.numeric(obs[2:3]), paramGP=paramGP, surface=S.vec, frame=gridFullL[1:(nrowcol^2),], 
                                nrowcol=nrowcol)
        # add nugget
        obs[4] <- rnorm(1, mean=obsT, sd=sqrt(sigmasq.nugget))
        path <- rbind(path,c(obs, tempGrad))
        gridFullL <- rbind(gridFullL, obs[2:3])
        S.vec <- c(S.vec, obsT)
        #print(path)
        betaF <- c(betaF, alpha[1])
      }
    }
  # this adds observation error (measurement error) if m>0. If m>0 then
  # no measurement noise
  path[,2] <- path[,2]+rnorm(nrow(path), 0, sqrt(m))
  path[,3] <- path[,3]+rnorm(nrow(path), 0, sqrt(m))
  path  <- cbind(path, betaF)
  #Include burn in period.
  if(burnIn>0){
    if(nrow(path)>burnIn){
      timeCut <- path[burnIn,1]
      path <- path[(burnIn:nrow(path)),]
      path[,1] <- path[,1] - timeCut
    }else{}
  }
  return(path) 
}
#################################################################################
# genPrefDataS ##################################################################
#################################################################################
# genPrefData creates numTracks number of tracks using parameters above
# output includes 4th column indicating current track number
# type=0 then movement based on temperature
# type=1 then movement based on temperature gradient (more realistic)
# start = 0 means uniform starting point (random), o/w multivariate normal with mean 0 and variance=start
# burnIn is number of burn in points (default 6). If burnIn=0 then no burn in
# timing is rate parameter for exponential distribution used to generate time difference between samples
genPrefDataHybridBehav <- function(n, movementParam, nrowcol, m, paramGP,
                                   numTracks, alpha=c(1,0), 
                                   rawDat=cbind(0,0,0),  start=0, burnIn=6,
                                   timing=10){
  param <- movementParam
  # Gaussian process parameters
  mean <- paramGP[1]
  phi <- paramGP[2]
  nugget <- paramGP[3]
  kappa <- paramGP[4]
  GPVar <- paramGP[5]
  nrow=nrowcol
  ncol=nrowcol
  x <- seq(-150, 150, length.out=nrow)
  y <- seq(-150, 150, length.out=ncol)
  nDim = nrow*ncol
  
  gridFull <- expand.grid(x,y)
  distMat <- as.matrix(dist(gridFull, method = "euclidean"))
  model <- RMwhittle(nu=kappa, var=GPVar, scale=phi)
  S <- matrix(slot(rawDat, 'data')$variable, ncol=1, nrow=(nrowcol^2))
  surface <- S + as.numeric(mean)
  # generate first track
  # make sure each track is of a certain length
  repeat{
    Dat <- genPrefTrackHybridBehav(n, param, paramGP, nrowcol, m, surface, gridFull, 
                                   model, rawDat, alpha, start, burnIn, timing)
    if (nrow(Dat) > (n/2)) break
  }    
  # label these tracks with a 1
  Dat <- cbind(Dat, 1)
  # if more than 1 track...
  if(numTracks>1){
    for(i in 2:numTracks){
      repeat{
        Dat1 <- genPrefTrackHybridBehav(n, param, paramGP, nrowcol, m, surface, gridFull, 
                                        model, rawDat, alpha,  start, burnIn, timing)
        # make sure each track is of a certain length
        if (nrow(Dat1) > (n/2)) break
      }
      #label each track
      Dat1 <- cbind(Dat1, i)
      Dat <- rbind(Dat,Dat1)
    }
  }else{}
  # convert data to a data frame
  Dat <- data.frame(Dat,row.names=NULL)
  return(list(Dat=Dat, surface=surface, x=x, y=y))
}
##########################################################################
# findTemp ###############################################################
##########################################################################
# this function finds the temperature at any given point within the bounds of the frame.
# either a weighted average of temperatures from the 4 vertices of the square around the point
# or 4 + n weighted temperatures for n observation locations within the square of interest
findTemp <- function(x, surface, frame, nrowcol){
  gridFull <- frame
  obs <- x
  xseq <- unique(gridFull[1:(nrowcol^2),1])
  yseq <- unique(gridFull[1:(nrowcol^2),2])
  # create rectangle/square around extrapolation point x
  tmp1 <- sign(obs[1] - xseq)
  tmp2 <- sign(obs[2] - yseq)
  n2x <- (1:(nrowcol-1))[(diff(tmp1)==(-2))] + 1
  n2y <- (1:(nrowcol-1))[(diff(tmp2)==(-2))] + 1
  n1x <- n2x - 1
  n1y <- n2y - 1
  p1 <- c(xseq[n1x], yseq[n1y])
  p2 <- c(xseq[n2x], yseq[n1y])
  p3 <- c(xseq[n1x], yseq[n2y])
  p4 <- c(xseq[n2x], yseq[n2y])
  # rectangle/square created as op
  op <- rbind(p1, p2, p3, p4)
  # see if any other observations in this area can be used to aid prediction
  gridSmall <- gridFull[(nrowcol^2):nrow(gridFull),]
  extraInit <- gridSmall[(gridSmall[,1]>p1[[1]] & gridSmall[,1]<p2[[1]] & gridSmall[,2]>p2[[2]] & gridSmall[,2]<p3[[2]]),]
  if(dim(extraInit)[1]>0){
    extra <- as.matrix(extraInit, ncol=2)
  }else{
    extra <- NULL
  }
  # add points for prediction to op
  op <- rbind(op, extra)
  # create a distance matrix of all vertices of square and observation in the middle
  finalSelec <- matrix(rbind(op, as.numeric(obs)), ncol=2)
  np <- dim(finalSelec)[1]
  distMat <- outer(as.numeric(obs)[1],as.numeric(finalSelec[,1]),"-")^2
  distMat <- distMat + outer(as.numeric(obs)[2],as.numeric(finalSelec[,2]),"-")^2
  distSquare <- sqrt(distMat)[1:(length(distMat)-1)]
  stdDistSquare <- (1/distSquare)/sum((1/distSquare))
  stdDistSquare[is.nan(stdDistSquare)] <- 1
  gridRowF <-  which(gridFull[,1] %in% as.numeric(op[,1]) & gridFull[,2] %in% op[,2])
  # return final predicted value
  obsT <- sum(surface[gridRowF] * stdDistSquare)
  return(obsT)
}
##########################################################################
# findTempGPV ############################################################
##########################################################################
# this function finds the temperature at any given point when generating data
# by 'filling' in the temeperature using conditional normal densities.
findTempGPV <- function(x, paramGP, surface, frame, nrowcol, neighbours){
  gridFull <- frame
  # location at which to find temperature
  newP <- x
  # spatial parameters  c(mean, phi, nugget, kappa, GPVar)
  mu <- paramGP[1] 
  # scale (range)
  phi <- abs(paramGP[2]) 
  # nugget variance
  sigmasq.nugget <- abs(paramGP[3]) 
  # marginal variance (partial sill)
  GPVar <- abs(paramGP[5]) 
  # kappa treated as known
  kappa <- paramGP[4]
  # make matrix of location
  xMat <- (matrix(newP, ncol=2))
  # find distance between location newP and all grid points
  dists <- as.numeric(sqrt((outer(xMat[,1],gridFull[,1],"-")^2) + (outer(xMat[,2],gridFull[,2],"-")^2)))
  # find closest distanced points
  closeDist <- sort(dists)[1:(neighbours)]
  # index these points
  indexNei <- which(dists%in%closeDist)
  # this is to make sure we only have 1 instance of newp in
  # the closest point list
  if(((nrowcol^2)+1)%in%indexNei){
    indexNei <- indexNei
  }else{
    indexNei <- c(which(dists%in%closeDist), ((nrowcol^2)+1))
  }
  # create list of all close locations and newP
  gridSmall <- rbind(gridFull[indexNei[-(neighbours+1)],], newP)
  # find all their respective distances
  distMatS <- sqrt((outer(gridSmall[,1],gridSmall[,1],"-")^2) + (outer(gridSmall[,2],gridSmall[,2],"-")^2))
  # find matern covariance matrix
  SigmaFullS <- geoR::matern(distMatS, phi, kappa)*GPVar
  newLength <- neighbours + 1
  # extract the chunks
  si12 <- SigmaFullS[1:neighbours, newLength]
  si11 <- SigmaFullS[1:neighbours, 1:neighbours]
  #si22 <- SigmaFullS[newLength, newLength]
  si21 <- t(si12)
  # extract the chunks
  mu <- rep(mu, newLength)
  mu1 <- mu[1:neighbours]
  mu2 <- mu[newLength]
  # generate a bivariate normal random vector
  # (X_4, X_5) ~ N( (mu[4], mu[5]), si22 )
  h1 <- rnorm(newLength)
  oldSurfaceDiff <- surface[indexNei[1:neighbours]] - mu1
  # now compute the mean of (X_1, X_2, X_3 | X_4, X_5)
  muc <- mu2 + si21 %*% solve(si11, oldSurfaceDiff)
  # new temperature prediction
  tempNew <- muc
  return(tempNew)
}

##########################################################################
# tempDrift###############################################################
##########################################################################
# find drift function/temp gradient at a particular location
# gp=1 if using findTempGp, gp=0 otherwise
tempDrift <- function(xy, paramGP, surface, frame, nrowcol, neighbours=50){
   mu <- -(jacobian(func=findTempGPV, x=xy, paramGP=paramGP, surface=surface, frame=frame, nrowcol=nrowcol,
                     neighbours=neighbours, method="simple")) 
 return(mu)
}