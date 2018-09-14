source('dataFncs.R')
compile("TMBfile.cpp")
dyn.load(dynlib("TMBfile"))
IGN <- function(pred, act, var) {
  ((pred - act)^2) / var + log(var)
}
# constant mean
mean <- 5
# scale (range)
phi <- 25 
# nugget (measurement) variance
nugget <- 0.1
# smoothness (assumed known in estimation)
kappa <- 2
# marginal variance (partial sill)
GPVar  <- 1.5
# define the covariance model
model <- RMwhittle(nu=kappa, var=GPVar, scale=phi)
# finally trend = 0 just means our mean trend is constant (mean from above)
trend <- 0

# is starting location random? (0 = yes and >0 multivariate normal with
# mean 0 and diagonal covariance matrix with variance=start)
start <- 0
# alpha[1] defines starting value of beta_1 from eq (3.5)
# alpha[2:3] are currently both equal to \alpha from eq (3.2). They could be changed to
# adopt preferential movement varying in strength across latitude/longitude.
alpha <- c(-1.5, 100, 100)
# the number of tracks in the simulation
numTracks <- 3
# the number of data points to simulate per track
n <- 360
# the number of observations to throw out per track (ie/ total sample
# size per track is n-burnIn)
burnIn <- 60
# the rate parameter of the exponential distribution used to generate the sampling times
timing <-  10
# measurement location noise (currently not included in models)
noise <- 0
# movement parameters
# behaviour (beta) standard deviation parameter (\sigma_{\beta} in eq (3.5))
behavSD <- .1 
# movement standard deviation parameter (diag(\Sigma) in eq (3.3))
moveSD <- 3 
# combine standard deviations for later use
dataParam <- c(behavSD, moveSD)


# define the domain to simulate data
# how man rows/columns
nrowcol <- 51
x <- seq(-150, 150, length.out=nrowcol)
y <- seq(-150, 150, length.out=nrowcol)
# simulation grid
gridFull <- expand.grid(x,y)

# l is the number of grid cells in each direction across our
# grid going from -150 to 150 degrees lat and lon for our model fitting
# and prediction.
l <- 26
xseq <- (seq(-150, 150, length.out=l))
yseq <- (seq(-150, 150, length.out=l))
# create the prediction lattice
lattice <- expand.grid(xseq,yseq)
colnames(lattice) <- c("Y1New", "Y2New")

# initiate objects
nonPrefParams <- NULL
prefParams <- NULL
postBias <- NULL
krigBias <- NULL
nonPrefBias <- NULL
postIGN <-  NULL
krigIGN <- NULL
nonPrefIGN <- NULL
nonPrefParams <- array(NA, dim=c(1, 4))
prefParams <- array(NA, dim=c(1, 8))
set.seed(6351)
# simulate the random field
rawDat <- RFsimulate(model, x=as.matrix(gridFull),  exactness=TRUE)
# simulate the observations and sampling locations
Xsim <- genPrefDataHybridBehav(n=n, movementParam=dataParam, nrowcol=nrowcol, m=0, 
                               paramGP=c(mean, phi, nugget, kappa, GPVar), numTracks=numTracks, 
                               alpha=alpha, rawDat=rawDat, start=start, burnIn = burnIn, timing=timing)
# extract sampling data
data <- Xsim$Dat
# extract true surface data
surface <- Xsim$surface
colnames(data) <- c("Time", "Lon", "Lat", "Temp", "gradientX", "gradientY", "Beta", "Track")
# here is how the data (locations and respective latent field measurements)
head(data)
# plot of generated data
image.plot(x,y,matrix(rawDat$variable1+mean, nrow=51, ncol=51),
           xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
totalPoints = 0
colourList <- c("purple", "blue", "green", "grey30", "black", "pink")
for(k in 1:numTracks){
  numPoints <- sum(Xsim$Dat$V8==k)
  points(Xsim$Dat$V2[(totalPoints+1):(totalPoints+numPoints)],
         Xsim$Dat$V3[(totalPoints+1):(totalPoints+numPoints)], type='b', pch=19, cex=.5, col=colourList[k])
  totalPoints <- totalPoints + numPoints
}
# now we thin the data to 300 locations in total for analysis 
selection <- seq(1, nrow(data), length.out = 300)
dataThin <- data[selection, ]
# plot of data to be analysed
image.plot(x,y,matrix(rawDat$variable1+mean, nrow=51, ncol=51),
           xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
totalPoints = 0
colourList <- c("purple", "blue", "green", "grey30", "black", "pink")
for(k in 1:numTracks){
  numPoints <- sum(dataThin$Track==k)
  points(dataThin$Lon[(totalPoints+1):(totalPoints+numPoints)],
         dataThin$Lat[(totalPoints+1):(totalPoints+numPoints)], type='b', pch=19, cex=.5, col=colourList[k])
  totalPoints <- totalPoints + numPoints
}
# replace data with thinned version
data=dataThin
# obtain sampling times
tsim <- data[,1]
# number of observations in total
numObs <- nrow(data)
# Generate random measurements
# create trackID which records when new tracks start in the dataframe
trackLength <- NULL
trackId <- 0
for(i in 1:numTracks){
  trackLength <- c(trackLength, length(which(data$Track==i)))
  trackId <- c(trackId, sum(trackLength))
}
# create a set of locations which allows for gradients to be calculated in cpp file
Yobs <- data$Temp
Y1New <- data[,2]
Y2New <- data[,3]
for(i in 1:length(data[,1])){
  Y1New <- c(Y1New,data[i,2]+.5, data[i,2])
  Y2New <- c(Y2New,data[i,3], data[i,3]+.5)
}
# combine prediction lattice with sampling locations and gradient locations
predGrid <-  rbind(cbind(Y1New, Y2New), lattice)

# create INLA mesh
mesh <- inla.mesh.create(loc = predGrid, extend = T, refine = T) 
# now create an index that matches sampling locations with mesh locations 
ii0 <- mesh$idx$loc
# Create data for TMB
dataTMB <- list(tsim=tsim,Y1=Y1New, Y2=Y2New, Y=Yobs, trackId=trackId,  meshidxloc=mesh$idx$loc-1)


# using SPDE method from R-INLA with alpha=2 (kappa=1)
dataTMB$spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
# create our own sparse precision matrix for alpha=3 (kappa=2)
M0 <- (inla.spde2.matern(mesh, alpha=2)$param.inla)["M0"]$M0
M1 <- (inla.spde2.matern(mesh, alpha=2)$param.inla)["M1"]$M1
M2 <- (inla.spde2.matern(mesh, alpha=2)$param.inla)["M2"]$M2
M3 <- M2%*%solve(M0)%*%M1
M3 <- as(M3, "dgTMatrix")
dataTMB$spde[["M3"]]<- M3
# number of rows in SPDE object
n_s = nrow(dataTMB$spde$M0)
# vector of 1's used in TMB (this should be updated)
dataTMB$Ind <- rep(1, n_s)
# create geodata object
obj1 <- cbind(cbind(dataTMB$Y1, dataTMB$Y2)[1:length(dataTMB$Y),], dataTMB$Y)
geodata <- as.geodata(obj1, coords.col = 1:2, data.col = 3)

## Parameter Estimation
standardMLE <- likfit(geodata, coords = geodata$coords, data = geodata$data, kappa=kappa, ini=c(.5,.5))

parameters <- list(
  S = rep(0, n_s),
  beta = rep(0, length(dataTMB$Y)),
  mu = standardMLE$beta,
  log_papertau = 3,
  log_kappa = log(1/standardMLE$phi),
  alpha = rnorm(1,alpha[2], 0.25),
  log_d = log(dataParam[2]),
  log_sdbehav = log(dataParam[1])
)
# create TMB object (note= random=c("S", "beta") to
# integrate out random field and latent behvaiour states)
obj <- MakeADFun(dataTMB, parameters, random=c("S", "beta"), DLL="TMBfile", method = "nlminb", hessian=FALSE, silent=T)
# conduct maximisation
opt <- try( nlminb(obj$par,obj$fn,obj$gr, control=list(rel.tol=1e-7)) )
# rerun up to 4 times in case of any gradient errors
for(m in 1:4){
  if(class(opt) != 'try-error' && opt$convergence == 0) {
    print("Success!")
  }
  else{ 
    paste0("Failed, try number ", m)
    lengthPar <- length(obj$env$last.par.best)
    tmp <- obj$env$last.par.best[(lengthPar-5):lengthPar] + 0.01
    opt <- try(nlminb(tmp,obj$fn,obj$gr, control=list(rel.tol=1e-7)))
  }
}
# Extract sigma^2 (partial sill)
report_spde <- obj$report()

# Obtain the standard errors    
sdre <- try( sdreport(obj) )
if( class(sdre) != 'try-error') {
  # input params
  summary(sdre, "fixed")
}
# prediction variance from TMB
predVar <- (summary(sdre, "random")[(length(Y1New)+1):(length(Y1New)+nrow(lattice)),2])^2

# conduct simple kriging using standard MLE plug-in parameters
SKDat <- krige.control(obj.model = standardMLE, type.krige = "SK")
# now predict at the prediction "lattice" locations where signal=T is used
# to specify that there was measurement error on our data observations
nonPredPref <- krige.conv(geodata, loc = lattice, krige = SKDat, output=list(signal=T))
# finally we obtain preferential model field prediction from mode of [S|Y,X]
modePred <- obj$env$last.par.best[(length(Y1New)+1):(length(Y1New)+nrow(lattice))]
# non-pref predictions
nonPrefPred <- nonPredPref$predict
# non-pref prediction variance
nonPrefVar <- nonPredPref$krige.var


# obtain real data on prediction lattice
# match indicies of full grid used to simulate data and prediction lattice
matchedIndic <- row.match(lattice,gridFull)
rawDatSmall <- rawDat$variable1[matchedIndic] + mean

IgnScorePost <- IGN(modePred, rawDatSmall, predVar)
IgnScoreNonPref <- IGN(nonPrefPred, rawDatSmall, nonPredPref$krige.var)
mean(IgnScorePost)
mean(IgnScoreNonPref)

# Before plots we will ignore locations far from observations (note that this is optional)
# set points far away from track as NA
sampLocs <- cbind(data$Lon, data$Lat)
colnames(sampLocs)=c("Y1New", "Y2New")
NAsampLocs <- rep(NA, nrow(lattice))
distMatLarge <- as.matrix(dist(rbind(lattice, sampLocs)))[(nrow(lattice)+1):(nrow(sampLocs) + nrow (lattice)), 1 : nrow(lattice)]
for(i in 1:nrow(lattice)){
  distMin <- min(as.matrix(dist(rbind(lattice[i,], sampLocs)))[1,-1])
  if(distMin<45){NAsampLocs[i] <- 1}
  else{}
}
setNA <- which(is.na(NAsampLocs))
modePredNA <- modePred
rawDatSmallNA <- rawDatSmall
predVarNA <- predVar
nonPrefPredNA <- nonPrefPred
nonPrefVarNA <- nonPrefVar
IgnScorePostNA <- IgnScorePost
IgnScoreNonPrefNA <- IgnScoreNonPref
modePredNA[setNA] <- NA
rawDatSmallNA[setNA] <- NA
predVarNA[setNA] <- NA
nonPrefPredNA[setNA] <- NA
nonPrefVarNA[setNA] <- NA
IgnScorePostNA[setNA] <- NA
IgnScoreNonPrefNA[setNA] <- NA

# Now we compare the true and predicted fields at non-NA locations
zlimPred=c(min(rawDatSmallNA,modePredNA,nonPrefPredNA, na.rm=T), max(rawDatSmallNA,modePredNA,nonPrefPredNA, na.rm=T))
image.plot(xseq,yseq,matrix(rawDatSmallNA, nrow=26, ncol=26),
           zlim=zlimPred, xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
totalPoints = 0
colourList <- c("purple", "blue", "green", "grey30", "black", "pink")
for(k in 1:numTracks){
  numPoints <- sum(Xsim$Dat$V8==k)
  points(Xsim$Dat$V2[(totalPoints+1):(totalPoints+numPoints)],
         Xsim$Dat$V3[(totalPoints+1):(totalPoints+numPoints)], type='b', pch=19, cex=.5, col=colourList[k])
  totalPoints <- totalPoints + numPoints
}
image.plot(xseq,yseq,matrix(modePredNA, nrow=26, ncol=26),
           zlim=zlimPred, xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
totalPoints = 0
colourList <- c("purple", "blue", "green", "grey30", "black", "pink")
for(k in 1:numTracks){
  numPoints <- sum(Xsim$Dat$V8==k)
  points(Xsim$Dat$V2[(totalPoints+1):(totalPoints+numPoints)],
         Xsim$Dat$V3[(totalPoints+1):(totalPoints+numPoints)], type='b', pch=19, cex=.5, col=colourList[k])
  totalPoints <- totalPoints + numPoints
}
image.plot(xseq,yseq,matrix(nonPrefPredNA, nrow=26, ncol=26),
           zlim=zlimPred, xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
totalPoints = 0
colourList <- c("purple", "blue", "green", "grey30", "black", "pink")
for(k in 1:numTracks){
  numPoints <- sum(Xsim$Dat$V8==k)
  points(Xsim$Dat$V2[(totalPoints+1):(totalPoints+numPoints)],
         Xsim$Dat$V3[(totalPoints+1):(totalPoints+numPoints)], type='b', pch=19, cex=.5, col=colourList[k])
  totalPoints <- totalPoints + numPoints
}
zlimIGN=c(min(IgnScorePostNA,IgnScoreNonPrefNA, na.rm=T), max(IgnScorePostNA,IgnScoreNonPrefNA, na.rm=T))
image.plot(xseq,yseq,matrix(IgnScorePostNA, nrow=26, ncol=26),
           zlim=zlimIGN, xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
totalPoints = 0
colourList <- c("purple", "blue", "green", "grey30", "black", "pink")
for(k in 1:numTracks){
  numPoints <- sum(Xsim$Dat$V8==k)
  points(Xsim$Dat$V2[(totalPoints+1):(totalPoints+numPoints)],
         Xsim$Dat$V3[(totalPoints+1):(totalPoints+numPoints)], type='b', pch=19, cex=.5, col=colourList[k])
  totalPoints <- totalPoints + numPoints
}
image.plot(xseq,yseq,matrix(IgnScoreNonPrefNA, nrow=26, ncol=26),
           zlim=zlimIGN, xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
totalPoints = 0
colourList <- c("purple", "blue", "green", "grey30", "black", "pink")
for(k in 1:numTracks){
  numPoints <- sum(Xsim$Dat$V8==k)
  points(Xsim$Dat$V2[(totalPoints+1):(totalPoints+numPoints)],
         Xsim$Dat$V3[(totalPoints+1):(totalPoints+numPoints)], type='b', pch=19, cex=.5, col=colourList[k])
  totalPoints <- totalPoints + numPoints
}
# Finally we can compare the IGN at the non-NA locations
mean(IgnScorePostNA, na.rm=T)
mean(IgnScoreNonPrefNA, na.rm=T)
