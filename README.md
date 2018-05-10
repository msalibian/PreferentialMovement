Preferential Sampling with moving monitoring stations
================

This repository is a companion resource to the paper "Modelling ocean temperatures from bio-probes under preferential sampling" (submitted), by Daniel Dinsdale and Matias Salibian-Barrera and it contains code to illustrate how to apply the methods discussed in that paper.

Introduction
------------

Below we illustrate how to model data obtained from sensor tags mounted on marine mammals which may have been preferentially obtained. The `R` code used below can be found in a single file here: [runfile.R](runfile.R). The required functions are in the file [dataFncs.R](dataFncs.R), and the negative log-likelihood function required by [TMB](https://cran.r-project.org/package=TMB) is in the file [TMBfile.cpp](TMBfile.cpp).

The example below follows that in Section 4 of the paper and uses the Preferential-CRW model to simulate preferentially sampled animal tracks and corresponding sea surface temperature observations. We then compare parameter estimation and prediction using standard methods and the preferential model in TMB.

### Warning

Note that running this code may require a large amount of RAM (we recommend 16GB).

Generating a data set
---------------------

First we source the file [dataFncs.R](dataFncs.R) which contains the necessary functions to generate the data.

``` r
source('dataFncs.R')
```

Now we specify the parameters to generate the data set. These can be altered to vary the properties of the latent field to be sampled and also to change the movement patterns of the sampler. Refer to the paper for more details on the model.

First we specify the field parameters sasuming a matern covariance structure:

``` r
# constant mean
mean <- 5
# scale (range)
phi <- 15 
# nugget (measurement) variance
nugget <- 0.1
# smoothness (assumed known in estimation)
kappa <- 2
# marginal variance (partial sill)
GPVar  <- 3
# define the covariance model
model <- RMwhittle(nu=kappa, var=GPVar, scale=phi)
# finally trend = 0 just means our mean trend is constant (mean from above)
trend <- 0
```

Next we specify the parameters that determine the movement/sampler properties:

``` r
# is starting location random? (0 = yes and >0 multivariate normal with
# mean 0 and diagonal covariance matrix with variance=start)
start <- 0
# alpha[1] defines starting value of beta_1 from eq (3.5)
# alpha[2:3] are currently both equal to \alpha from eq (3.2). They could be changed to
# adopt preferential movement varying in strength across latitude/longitude.
alpha <- c(.5, 150, 150) 
# the number of tracks in the simulation
numTracks <- 4
# the number of data points to simulate per track
n <- 80
# the number of observations to throw out per track (ie/ total sample
# size per track is n-burnIn)
burnIn <- 20
# measurement location noise (currently not included in models)
noise <- 0
# movement parameters
# behaviour (beta) standard deviation parameter (\sigma_{\beta} in eq (3.5))
behavSD <- .1 
# movement standard deviation parameter (diag(\Sigma) in eq (3.3))
moveSD <- 12 
# combine standard deviations for later use
dataParam <- c(behavSD, moveSD)
```

We now create a lattice for data simulation and also for model fitting/predictions. These can be the same if `nrowcol` = `l` but we can choose different grid sizes for computational efficiency.

``` r
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
```

<!-- ```{r initiate, include=FALSE} -->
<!-- # initiate objects -->
<!-- # nonPrefParams <- NULL -->
<!-- # prefParams <- NULL -->
<!-- # postBias <- NULL -->
<!-- # krigBias <- NULL -->
<!-- # nonPrefBias <- NULL -->
<!-- # postIGN <-  NULL -->
<!-- # krigIGN <- NULL -->
<!-- # nonPrefIGN <- NULL -->
<!-- # nonPrefParams <- array(NA, dim=c(1, 4)) -->
<!-- # prefParams <- array(NA, dim=c(1, 8)) -->
<!-- ``` -->
Now we can generate the data. We first simulate the latent field and then run the sampler using `genPrefDataHybridBehav` which can be found in `dataFncs.R`. From here we will extract the data and the so-called true surface on the lattice and observed locations

``` r
# simulate the random field
set.seed(1)
rawDat <- RFsimulate(model, x=as.matrix(gridFull),  exactness=TRUE)
# simulate the observations and sampling locations
Xsim <- genPrefDataHybridBehav(n=n, movementParam=dataParam, nrowcol=nrowcol, m=0, 
                               paramGP=c(mean, phi, nugget, kappa, GPVar), numTracks=numTracks, 
                               alpha=alpha, rawDat=rawDat, start=start, burnIn = burnIn)
# extract sampling data
data <- Xsim$Dat
# extract true surface data
surface <- Xsim$surface
colnames(data) <- c("Time", "Lon", "Lat", "Temp", "Beta", "Track")
# here is how the data (locations and respective latent field measurements)
head(data)
```

    ##        Time        Lon        Lat     Temp      Beta Track
    ## 1 0.0000000  -4.109334  -85.63149 5.085735 0.5073743     1
    ## 2 0.2139963 -30.402035 -109.75711 6.749778 0.5005834     1
    ## 3 0.2189868 -29.278130 -108.76321 6.382570 0.5008400     1
    ## 4 0.2661443 -25.649978 -105.95862 6.503200 0.4983315     1
    ## 5 0.3997398 -10.621576  -95.71326 4.955609 0.5221208     1
    ## 6 0.5041566  -4.527735  -87.21801 4.935147 0.5091465     1

Here is how the data looks. Each colour is a different track and dots are sampling locations which are superimposed onto the unknown latent field.

![](README_files/figure-markdown_github/plotdata-1.png)

Now we compile the file [TMBfile.cpp](TMBfile.cpp) to use in TMB. This file contains the negative joint log-likelihood function −log(\[*X*, *Y*, *S*\]). Note that you must have installed the [TMB](https://cran.r-project.org/package=TMB) `R` package from CRAN.

``` r
compile("TMBfile.cpp")
dyn.load(dynlib("TMBfile"))
```

Next is some house keeping to prepare the data for TMB

``` r
# obtain sampling times
tsim <- data[,1]
# number of observations in total
numObs <- nrow(data)
# Generate random measurements
# create trackID which records when new tracks start in the dataframe
trackLength <- NULL
trackId <- 0
for(i in 1:numTracks){
  trackLength <- c(trackLength, length(which(data[,6]==i)))
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
```

Next we create a mesh using `inla.mesh.create` for the SPDE approach of `R-INLA`. We mush be careful to specify an index that matches sampling locations with mesh locations, but also change indexing for use in `C++`.

``` r
# create INLA mesh
mesh <- inla.mesh.create(loc = predGrid, extend = T, refine = T)
# now create an index that matches sampling locations with mesh locations
ii0 <- mesh$idx$loc
# Create data for TMB
dataTMB <- list(tsim=tsim,Y1=Y1New, Y2=Y2New, Y=Yobs, trackId=trackId,  meshidxloc=mesh$idx$loc-1)
```

Now we will create our sparse precision matrix for smoothness (*κ*) 2, which enables the field to be differentiable in the mean square sense. For details on this part see Appendix A.

``` r
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
```

Parameter Estimation
--------------------

Time to fit some models! First let us fit a standard model using `likfit` from the `geoR` package. This ignores any preferential effect and conditions on the sampling locations *X*.

``` r
standardMLE <- likfit(geodata, coords = geodata$coords, data = geodata$data, kappa=kappa, ini=c(.5,.5))
```

    ## ---------------------------------------------------------------
    ## likfit: likelihood maximisation using the function optim.
    ## likfit: Use control() to pass additional
    ##          arguments for the maximisation function.
    ##         For further details see documentation for optim.
    ## likfit: It is highly advisable to run this function several
    ##         times with different initial values for the parameters.
    ## likfit: WARNING: This step can be time demanding!
    ## ---------------------------------------------------------------
    ## likfit: end of numerical maximisation.

``` r
(standardMLE)
```

    ## likfit: estimated model parameters:
    ##      beta     tausq   sigmasq       phi 
    ## " 4.4319" " 0.0739" " 2.7101" "15.8464" 
    ## Practical Range with cor=0.05 for asymptotic range: 85.06966
    ## 
    ## likfit: maximised log-likelihood = -122.8

Next we will fit the model in `TMB`. First we define the parameters for the model (including latent states). Our latent states are the field *S* and behavioural states *b**e**t**a*'s. We start the other parameters

``` r
parameters <- list(
  S = rep(0, n_s),
  beta = rep(0.5, length(dataTMB$Y)),
  mu = standardMLE$beta,
  log_papertau = 3,
  log_kappa = log(1/standardMLE$phi),
  log_tau = log(standardMLE$tausq),
  alpha = rnorm(1,alpha[2], 0.25),
  log_d = log(dataParam[2]),
  log_sdbehav = log(dataParam[1])
)
# create TMB object (note= random=c("S", "beta") to
# integrate out random field and latent behvaiour states)
obj <- MakeADFun(dataTMB, parameters, random=c("S", "beta"), DLL="TMBfile", method = "nlminb", hessian=FALSE, silent=T)
# conduct maximisation
opt <- try( nlminb(obj$par,obj$fn,obj$gr, control=list(rel.tol=1e-7)) )
# Extract sigma^2 (partial sill)
report_spde <- obj$report()
```

``` r
# check convergence
opt$convergence
```

    ## [1] 0

``` r
# Obtain the standard errors
sdre <- try( sdreport(obj) )
if( class(sdre) != 'try-error') {
  # input params
  summary(sdre, "fixed")
}
```

    ##                Estimate  Std. Error
    ## mu             4.939458  0.70357739
    ## log_papertau   3.539569  0.16696319
    ## log_kappa     -2.869189  0.14946801
    ## log_tau       -1.285581  0.05806115
    ## alpha        174.620945 23.13155555
    ## log_d          2.574603  0.03588378
    ## log_sdbehav   -1.870642  0.31034349

``` r
# prediction variance from TMB
predVar <- (summary(sdre, "random")[(length(Y1New)+1):(length(Y1New)+nrow(lattice)),2])^2
```

Prediction
----------

Now we have obtained parameter estimates for the standard method and for the preferential model using `TMB`. To predict using the non-preferential model we will use kriging with plug-in parameters obtained from the standard `likfit` function. For the preferential model we use the mode of the \[*S*|*Y*, *X*\] at the optimal parameters. This is provided by *T**M**B* as part of the Laplace approximation procedure and is defined in eq (2.7).

``` r
# conduct simple kriging using standard MLE plug-in parameters
SKDat <- krige.control(obj.model = standardMLE, type.krige = "SK")
# now predict at the prediction "lattice" locations where signal=T is used
# to specify that there was measurement error on our data observations
nonPredPref <- krige.conv(geodata, loc = lattice, krige = SKDat, output=list(signal=T))
```

    ## krige.conv: model with constant mean
    ## krige.conv: Kriging performed using global neighbourhood

``` r
# finally we obtain preferential model field prediction from mode of [S|Y,X]
modePred <- obj$env$last.par.best[(length(Y1New)+1):(length(Y1New)+nrow(lattice))]
# non-pref predictions
nonPrefPred <- nonPredPref$predict
# non-pref prediction variance
nonPrefVar <- nonPredPref$krige.var
```

Next we want to be able to compare these predictions to the real values of the field at the prediction points.

``` r
# obtain real data on prediction lattice
# match indicies of full grid used to simulate data and prediction lattice
matchedIndic <- row.match(lattice,gridFull)
rawDatSmall <- rawDat$variable1[matchedIndic] + mean
```

Now let us calculate the mean ignorance score for each method on this data set (MIGN from eq (4.2)). Recall that the ignorance function (IGN) is given by

``` r
IGN <- function(pred, act, var) {
  ((pred - act)^2) / var + log(var)
}
```

Then the MIGN can be computed as follows:

``` r
IgnScorePost <- IGN(modePred, rawDatSmall, predVar)
IgnScoreNonPref <- IGN(nonPredPref$predict, rawDatSmall, nonPredPref$krige.var)
mean(IgnScorePost)
```

    ## [1] 1.170347

``` r
mean(IgnScoreNonPref)
```

    ## [1] 1.299883

Finally we can plot the IGN scores and compare predictive surfaces from the non-preferential and preferential models. We consider only prediction locations in regions near the sampling locations: ![](README_files/figure-markdown_github/showign-1.png)![](README_files/figure-markdown_github/showign-2.png)![](README_files/figure-markdown_github/showign-3.png) Note that the mean IGN for the following two plots are 0.41 (TMB) and 0.64 (kriging) respectively. ![](README_files/figure-markdown_github/plotign2-1.png)![](README_files/figure-markdown_github/plotign2-2.png)
