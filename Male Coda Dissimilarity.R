# Male codas
library(spdep)

## Calculate geographic distance

male.gps.data <- read.csv('male.coda.gps.csv')
gps.data.updated <- merge(combined.codas.all.sites,male.gps.data, by.x = 'individual', by.y = 'pair')
gps.data.updated <- droplevels(subset(gps.data.updated, site=="SA"))

gps.data.updated[which(is.na(gps.data.updated))]


geo_distance <- spDists(cbind(gps.data.updated$lon,gps.data.updated$lat),longlat = T)
geo_matrix <- as.matrix(geo_distance)
(geo_matrix[1,])

geo_matrix2 <- geo_matrix
geo_matrix2[lower.tri(geo_matrix2 , diag=TRUE)] <- NA
geo_vec <- na.exclude(as.vector(geo_matrix2))


##Calculate acoustic distance
#fit.pca <- princomp(features[6:28], scores=T)#
pca.dists <- dist(gps.data.updated[,3:24])
str(pca.dists)

acoustic.dist <- pca.dists
acou_matrix <- as.matrix(acoustic.dist)

aco_matrix2 <- acou_matrix
aco_matrix2[lower.tri(aco_matrix2, diag=TRUE)] <- NA
aco_vec <- na.exclude(as.vector(aco_matrix2))

residual_vectors <- cbind(aco_vec, geo_vec)
resid_vecs <- as.data.frame(residual_vectors)
str(resid_vecs)

resid_vecs <- subset(resid_vecs, subset=geo_vec >0)


### Calculate confidence intervals 
acoustic.dist.dataframe<- cbind.data.frame(resid_vecs$geo_vec,resid_vecs$aco_vec)

new.acoustic.data.frame <- subset(acoustic.dist.dataframe, resid_vecs$geo_vec>=0.5)

x <- new.acoustic.data.frame$`resid_vecs$geo_vec`
y <- new.acoustic.data.frame$`resid_vecs$aco_vec`


my.spar <- 1.5 ## this is the smoothing parameter for all spline bootstraps
sp.frame <- data.frame(x=x,y=y)
sp.resampler <- function() {
  n <- nrow(sp.frame)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(sp.frame[resample.rows,])
}

sp.spline.estimator <- function(data,m=300) {
  # Fit spline to data, with cross-validation to pick lambda
  fit <- smooth.spline(x=data[,1],y=data[,2],spar=my.spar)
  # Set up a grid of m evenly-spaced points on which to evaluate the spline
  eval.grid <- seq(from=min(x),to=max(x),length.out=m)
  # Slightly inefficient to re-define the same grid every time we call this,
  # but not a big overhead
  # Do the prediction and return the predicted values
  return(predict(fit,x=eval.grid)$y) # We only want the predicted values
}

sp.spline.cis <- function(B,alpha,m=300) {
  spline.main <- sp.spline.estimator(sp.frame,m=m)
  # Draw B boottrap samples, fit the spline to each
  spline.boots <- replicate(B,sp.spline.estimator(sp.resampler(),m=m))
  # Result has m rows and B columns
  cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2)
  cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
              x=seq(from=min(x),to=max(x),length.out=m)))
}

sp.cis <- sp.spline.cis(B=1000,alpha=0.05)

plot(x,y,xlab="Geographic Distance (km)",
     ylab="Male Coda Dissimilarity",col= "grey80", pch =".", cex=.7)
#smooth.spline(y,x, cv=FALSE)
lines(x=sp.cis$x,y=sp.cis$main.curve, col="red")
lines(x=sp.cis$x,y=sp.cis$lower.ci)
lines(x=sp.cis$x,y=sp.cis$upper.ci)

