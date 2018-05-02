# process 2 small-simulation outputs
# 
library(fields)
library(ggplot2)
# Load data
load('githubtest.RData')

tmp1 <- complete.cases(nonPrefParams) & complete.cases(prefParams)
tmp1 <- tmp1 & (prefParams[,8] == 0)
#tmp1[155] = FALSE
# pref
# tmp1[c(9,31, 57, 80 )] = FALSE
# non-pref
# tmp1[which(is.na(IgnScorePost0))[1:100]] = FALSE



npr <- apply(IgnScoreNonPref0[tmp1,,], 1, mean)#[1:100]
pr <- apply(IgnScorePost0[tmp1,,], 1, mean)#[1:100]
prefParams1=prefParams[tmp1,]#[1:100,]
nonPrefParams1=nonPrefParams[tmp1,]#[1:100,]
TrueSurface1=TrueSurface[tmp1,,]#[1:100,,]
PredPost01=PredPost0[tmp1,,]#[1:100,,]
PredNonPref01=PredNonPref0[tmp1,,]#[1:100,,]

# Plot parameters IGN scores
##########################################################
# IGN
boxplot(pr, npr, names=c('Pref', 'NonPref'))
# IGN difference (gives clearer picture of performance I think?)
boxplot(pr - npr)
abline(h=0)
# Parameters
boxplot(prefParams1[,1], nonPrefParams1[, 1])
abline(h=5)
boxplot(prefParams1[,2], nonPrefParams1[, 2])
abline(h=15)
boxplot(prefParams1[,3], nonPrefParams1[, 3])
abline(h=3)
boxplot(prefParams1[,4], nonPrefParams1[, 4])
abline(h=0.1)
# Plot average RMSPE of each lcoation across all sims
# (currently sets negative values to 0)
prMSPE <- sqrt((TrueSurface1 - PredPost01)^2)
nprMSPE <- sqrt((TrueSurface1 - PredNonPref01)^2)
diffMSPE <- abs(nprMSPE)-abs(prMSPE)
avgDiffMSPE <- apply(diffMSPE, c(2,3), mean)
badPred <- which(avgDiffMSPE<0)
avgDiffMSPE[badPred] <- NA
# Do plot
x <- seq(-150, 150, length.out=26)
y <- seq(-150, 150, length.out=26)
image.plot(x,y,avgDiffMSPE, xlab="Longitude", ylab="Latitude", col=heat.colors(10))
# Plot IGN map
IGNdiffMap <- apply(IgnScorePost0[tmp1,,]-IgnScoreNonPref0[tmp1,,], c(2,3), mean)
badPredIGN <- which(IGNdiffMap>0)
IGNdiffMap[badPredIGN] <- NA
image.plot(x,y,IGNdiffMap, xlab="Longitude", ylab="Latitude", col=rev(heat.colors(10)))
#######################################################################################################################################
# Plots for paper #####################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Mean parameter ######################################################################################################################
#######################################################################################################################################
n <- nrow(nonPrefParams1)

#box plots of mean parameter estimates
paramLong <- matrix(NA, nrow=(4*n), ncol=3)
# paramLong[1:n, 1] <- prefParamsINLA[,1]
# paramLong[1:n, 2] <- "mu"
# paramLong[1:n, 3] <- "INLA"

paramLong[(n+1):(2*n), 1] <- prefParams1[,1]
paramLong[(n+1):(2*n), 2] <- "mu"
paramLong[(n+1):(2*n), 3] <- "Pref"

paramLong[(2*n+1):(3*n), 1] <- nonPrefParams1[,1]
paramLong[(2*n+1):(3*n), 2] <- "mu"
paramLong[(2*n+1):(3*n), 3] <- "NonPref"
# 
# paramLong[(3*n+1):(4*n), 1] <- prefParamsDig[,1]
# paramLong[(3*n+1):(4*n), 2] <- "mu"
# paramLong[(3*n+1):(4*n), 3] <- "MC2010"

paramLong <- data.frame(paramLong)

paramLong[,1] <- as.numeric(as.character(paramLong[,1]))

colnames(paramLong) <- c("value", "param", "method")

# pdf(file='C:/PDFs/meanParamEstimatesS2sim100.pdf')
p1 <- ggplot(na.omit(paramLong), aes(factor(method), value))
p1 + geom_boxplot(aes(fill = method)) + 
  geom_hline(aes(yintercept=5), colour="blue", alpha=.5, size = 2) +
  # ggtitle("Scale Parameter Estimates") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Parameter Estimate") + theme(legend.position="none", axis.text.y = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=30))
# dev.off()
#######################################################################################################################################
# Scale parameter #####################################################################################################################
#######################################################################################################################################
n <- nrow(nonPrefParams1)
paramLong <- matrix(NA, nrow=(13*n), ncol=4)

paramLong[(4*n+1):(5*n), 1] <- prefParams1[,2]
paramLong[(4*n+1):(5*n), 2] <- "phi"
paramLong[(4*n+1):(5*n), 4] <- 15

paramLong[(5*n+1):(6*n), 1] <- prefParams1[,3]
paramLong[(5*n+1):(6*n), 2] <- "sigma[2]"
paramLong[(5*n+1):(6*n), 4] <- 3


paramLong[(6*n+1):(7*n), 1] <- prefParams1[,4]
paramLong[(6*n+1):(7*n), 2] <- "nugget"
paramLong[(6*n+1):(7*n), 4] <- 0.1


paramLong[(7*n+1):(8*n), 1] <- nonPrefParams1[,2]
paramLong[(7*n+1):(8*n), 2] <- "phi"
paramLong[(7*n+1):(8*n), 4] <- 15

paramLong[(8*n+1):(9*n), 1] <- nonPrefParams1[,3]
paramLong[(8*n+1):(9*n), 2] <- "sigma[2]"
paramLong[(8*n+1):(9*n), 4] <- 3


paramLong[(9*n+1):(10*n), 1] <- nonPrefParams1[,4]
paramLong[(9*n+1):(10*n), 2] <- "nugget"
paramLong[(9*n+1):(10*n), 4] <- 0.1


paramLong[(4*n+1):(7*n),3] <- "Pref"
paramLong[(7*n+1):(10*n),3] <- "NonPref"


paramLong <- data.frame(paramLong)

paramLong[,1] <- as.numeric(as.character(paramLong[,1]))
paramLong[,4] <- as.numeric(as.character(paramLong[,4]))

colnames(paramLong) <- c("value", "param", "method", "true")


# plot phi
# pdf(file='C:/PDFs/phiEstimatesS2sim100.pdf')
p2 <- ggplot(na.omit(subset(paramLong, (param == "phi"))), aes(factor(method), value))
p2 + geom_boxplot(aes(fill = method)) + 
  geom_hline(aes(yintercept=true), colour="blue", alpha=.5, size = 2) +
  # ggtitle("Scale Parameter Estimates") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Parameter Estimate") + theme(legend.position="none", axis.text.y = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=30))
# dev.off()
# plot var
# pdf(file='C:/PDFs/varEstimatesS2sim100.pdf')
p3 <- ggplot(na.omit(subset(paramLong, (param == "sigma[2]"))), aes(factor(method), value))
p3 + geom_boxplot(aes(fill = method)) + 
  geom_hline(aes(yintercept=true), colour="blue", alpha=.5, size = 2) +
  # ggtitle("Scale Parameter Estimates") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Parameter Estimate") + theme(legend.position="none", axis.text.y = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=30))
# dev.off()
# plot nugget
# pdf(file='C:/PDFs/nuggetEstimatesS2sim100.pdf')
p4 <- ggplot(na.omit(subset(paramLong, (param == "nugget"))), aes(factor(method), value))
p4 + geom_boxplot(aes(fill = method)) + 
  geom_hline(aes(yintercept=true), colour="blue", alpha=.5, size = 2) +
  # ggtitle("Scale Parameter Estimates") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Parameter Estimate") + theme(legend.position="none", axis.text.y = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=30))
# dev.off()
#######################################################################################################################################
#######################################################################################################################################
# Movement parameters #################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
n <- nrow(nonPrefParams1)
paramLong <- matrix(NA, nrow=(4*n), ncol=4)

paramLong[(n+1):(2*n), 1] <- prefParams1[,5]
paramLong[(n+1):(2*n), 2] <- "alpha"
paramLong[(n+1):(2*n), 4] <- 150

paramLong[(2*n+1):(3*n), 1] <- prefParams1[,6]
paramLong[(2*n+1):(3*n), 2] <- "diffusion"
paramLong[(2*n+1):(3*n), 4] <- 12

paramLong[(3*n+1):(4*n), 1] <- prefParams1[,7]
paramLong[(3*n+1):(4*n), 2] <- "betaVar"
paramLong[(3*n+1):(4*n), 4] <- 0.1



paramLong[1:(4*n),3] <- "Pref"



paramLong <- data.frame(paramLong)

paramLong[,1] <- as.numeric(as.character(paramLong[,1]))
paramLong[,4] <- as.numeric(as.character(paramLong[,4]))

colnames(paramLong) <- c("value", "param", "method", "true")

# pdf(file='C:/PDFs/alphaParamEstimatesS2sim100.pdf')
p5 <- ggplot(na.omit(subset(paramLong, (param == "alpha"))), aes(factor(method), value))
p5 + geom_boxplot(aes(fill = method)) + 
  geom_hline(aes(yintercept=true), colour="blue", alpha=.5, size = 2) +
  # ggtitle("Scale Parameter Estimates") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Parameter Estimate") + theme(legend.position="none", axis.text.y = element_text(face="bold", size=15))
# dev.off()
# pdf(file='C:/PDFs/diffusionParamEstimatesS2sim100.pdf')
p6 <- ggplot(na.omit(subset(paramLong, (param == "diffusion"))), aes(factor(method), value))
p6 + geom_boxplot(aes(fill = method)) + 
  geom_hline(aes(yintercept=true), colour="blue", alpha=.5, size = 2) +
  # ggtitle("Scale Parameter Estimates") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Parameter Estimate") + theme(legend.position="none", axis.text.y = element_text(face="bold", size=15))
# dev.off()
# pdf(file='C:/PDFs/betaVarParamEstimatesS2sim100.pdf')
p7 <- ggplot(na.omit(subset(paramLong, (param == "betaVar"))), aes(factor(method), value))
p7 + geom_boxplot(aes(fill = method)) + 
  geom_hline(aes(yintercept=true), colour="blue", alpha=.5, size = 2) +
  # ggtitle("Scale Parameter Estimates") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Parameter Estimate") + theme(legend.position="none", axis.text.y = element_text(face="bold", size=15))
# dev.off()
#######################################################################################################################################
#######################################################################################################################################
# Ignorance Scores ####################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
n <- nrow(nonPrefParams1)
paramLong <- matrix(NA, nrow=(3*n), ncol=2)

paramLong[(n+1):(2*n), 1] <- pr - npr
paramLong[(n+1):(2*n), 2] <- "IGN Difference"


paramLong <- data.frame(paramLong)

paramLong[,1] <- as.numeric(as.character(paramLong[,1]))

colnames(paramLong) <- c("value", "method")
# pdf(file='C:/PDFs/IGNS2sim100.pdf')
p7 <- ggplot(na.omit(paramLong), aes(factor(method), value))
p7 + geom_boxplot(aes(fill = method)) + 
  # ggtitle("Scale Parameter Estimates") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("Mean Ignorance Score") + theme(legend.position="none", axis.text.y = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=30))
# dev.off()
#######################################################################################################################################
#######################################################################################################################################
# RMSPE Difference ####################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
prMSPE <- sqrt((TrueSurface1 - PredPost01)^2)
nprMSPE <- sqrt((TrueSurface1 - PredNonPref01)^2)
diffMSPE <- abs(prMSPE)-abs(nprMSPE)
avgDiffMSPE <- apply(diffMSPE, c(2,3), mean)
# badPred <- which(avgDiffMSPE<0)
# avgDiffMSPE[badPred] <- NA
# pdf(file='C:/PDFs/MSPEComparisonS2sim100.pdf')
#plot.new()
x <- seq(-150, 150, length.out=26)
y <- seq(-150, 150, length.out=26)
colfunc <- colorRampPalette(c("blue", "white", "red"))
myBreaks <- c(seq(min(avgDiffMSPE), 0, length.out=ceiling(1000/2) + 1), 
              -seq(-0.0003, min(avgDiffMSPE), length.out=ceiling(1000/2)))
image.plot(x,y,matrix(avgDiffMSPE, nrow=length(x), ncol=length(x)), 
           breaks=myBreaks, xlab="Longitude", ylab="Latitude", col=colfunc(1000))
# dev.off()
# # Plot IGN map
IGNdiffMap <- apply(IgnScorePost0[tmp1,,][1:100,,]-IgnScoreNonPref0[tmp1,,][1:100,,], c(2,3), mean)
badPredIGN <- which(IGNdiffMap>0)
# IGNdiffMap[badPredIGN] <- NA
# pdf(file='C:/PDFs/IGNComparisonS2sim100.pdf')
myBreaks <- c(seq(min(IGNdiffMap), 0, length.out=ceiling(1000/2) + 1), 
              -seq(-0.0003, min(IGNdiffMap), length.out=ceiling(1000/2)))
image.plot(x,y,IGNdiffMap, breaks=myBreaks, 
           xlab="Longitude", ylab="Latitude", col=colfunc(1000))
# dev.off()
