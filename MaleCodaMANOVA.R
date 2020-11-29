### R Code to investigate geographic variation

### MANOVA Model Code 
### Load required libraries 
library(corrplot)
library(bayesplot)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())
library(stringr)

#system("killall R")

## Load data 
#combined.codas.all.sites <- malecodadf#read.csv("combined.codas.all.sites.csv")
colnames(combined.codas.all.sites)
table(combined.codas.all.sites$individual)
length(unique(combined.codas.all.sites$individual))

ggpubr::ggdensity(data=combined.codas.all.sites,x="noterate",fill = "site")+
  scale_fill_manual(values=matlab::jet.colors(length(unique(combined.codas.all.sites$site))))

ggpubr::ggdensity(data=combined.codas.all.sites,x="meanbw",fill = "site")+
  scale_fill_manual(values=matlab::jet.colors(length(unique(combined.codas.all.sites$site))))

ggpubr::ggboxplot(data=combined.codas.all.sites,x='individual',y='range.bw')

hist(combined.codas.all.sites$call.dur)
table(combined.codas.all.sites$site)
hist((combined.codas.all.sites$minbw))

## Isolate relevant features from data set
# maxbw and max95 are highly correlated (0.9) so if we know max95 we have a pretty good idea of maxbw
# rest dur and call.dur also correlated (0.8)
d.manova <- combined.codas.all.sites[,c("rest.dur", "minbw","max95","noterate")]

## Check the structure of the data
cor(d.manova)

## Log transform data
d.manova <- log(d.manova)

pairs(d.manova)

## Check the structure of the data
str(d.manova)

### Set-up data to pass to Stan. 
# Integer-coded vector of group IDs
group.int <- as.numeric(combined.codas.all.sites$individual)

# Check structure 
table(group.int)

# Integer-coded vector of site IDs
site.int <- as.numeric(as.factor(combined.codas.all.sites$site))

# Check structure 
table(site.int)

# Center data matrix at feature means
col.means <- apply(d.manova, MARGIN=2, FUN="mean")
y.centered <- sweep(d.manova, MARGIN=2, STATS=col.means)

# Create a data list to pass to Stan
data_list <- list(
  K = dim(d.manova)[2],
  J= length(unique(combined.codas.all.sites$individual)),
  M= length(unique(combined.codas.all.sites$site)),
  N= dim(d.manova)[1],
  y= as.matrix(y.centered), ## features centered at zero
  group= group.int,
  site= site.int
)


# Code to run the STAN model
# NOTE: the .stan file must be linked in the code below
mfinal.stan = stan(file="MANOVA.male codas.stan", 
                   model_name = "mfinal", 
                   data=data_list, iter=3000, warmup=1500, chains=2, 
                   cores=2, 
                   control = list(stepsize = 0.5, adapt_delta = 0.99, max_treedepth = 20))

# Optional code to save the output
#save(mfinal.stan, file = "/Volumes/DJC HardDrive/stan.model.output.oct.2017.rda")

## Check model output
# Create traceplots to check for mixing for site level variance
draws <- as.array(mfinal.stan, pars="ICC_site")
mcmc_trace(draws)
stan_dens(mfinal.stan, pars=c("ICC_site"))
round(summary(mfinal.stan, pars=c("ICC_site"))$summary, 3)

# Create traceplots to check for mixing for group level variance
draws <- as.array(mfinal.stan, pars="ICC_group")
mcmc_trace(draws)
stan_dens(mfinal.stan, pars=c("ICC_group"))
round(summary(mfinal.stan, pars=c("ICC_group"))$summary, 3)

# Check degrees of freedom parameter
stan_dens(mfinal.stan, pars=c("DF_obs"))
round(summary(mfinal.stan, pars=c("DF_obs"))$summary, 3)


stan_dens(mfinal.stan, pars=c("DF_group"))
round(summary(mfinal.stan, pars=c("DF_group"))$summary, 3)


stan_dens(mfinal.stan, pars=c("DF_site"))
round(summary(mfinal.stan, pars=c("DF_site"))$summary, 3)

## Extract site-specific random intercepts for max.freq rate
site.intercepts <- extract(mfinal.stan, pars="site_rand_intercept")$site_rand_intercept
max.freq.intercepts <- data.frame(site.intercepts[, , 1]) # Intercepts for the 8th feature. 
names(max.freq.intercepts) <-levels(as.factor(combined.codas.all.sites$site))
str(max.freq.intercepts)

## Convert data into dataframe to pass to ggplot
CR <- cbind.data.frame(max.freq.intercepts$CR, rep("Crocker (CR)",length(max.freq.intercepts$CR)))
colnames(CR) <- c("samples","site")
DK <- cbind.data.frame(max.freq.intercepts$DK, rep("Deramakot (DK)",length(max.freq.intercepts$DK)))
colnames(DK) <- c("samples","site")
DV <- cbind.data.frame(max.freq.intercepts$DV, rep("Danum (DV)",length(max.freq.intercepts$DV)))
colnames(DV) <- c("samples","site")
IC <- cbind.data.frame(max.freq.intercepts$IC, rep("Imbak (IK)",length(max.freq.intercepts$IC)))
colnames(IC) <- c("samples","site")
KB <- cbind.data.frame(max.freq.intercepts$KB, rep("Kinabatangan (KB)",length(max.freq.intercepts$KB)))
colnames(KB) <- c("samples","site")
MB <- cbind.data.frame(max.freq.intercepts$MB, rep("Maliau (MB)",length(max.freq.intercepts$MB)))
colnames(MB) <- c("samples","site")
SAF <- cbind.data.frame(max.freq.intercepts$SA, rep("Kalabakan (KL)",length(max.freq.intercepts$SA)))
colnames(SAF) <- c("samples","site")
random.intercept.df.maxfreq <- rbind.data.frame(CR,DK, DV, IC, KB,MB, SAF)

# Check structure of data frame
head(random.intercept.df.maxfreq)
dim(random.intercept.df.maxfreq)

## Create color palette for figure 
my.cols <- matlab::jet.colors(7)


## Create Figure 5 in ggplot
ggplot(random.intercept.df.maxfreq, aes(x=samples, fill=site))+ geom_density(alpha=.8, bw=.015)+
  xlab("Deviation from Base Line")+ylab("")+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                                  axis.text.x  = element_text(size=20))+
  #scale_fill_brewer(palette="OrRd")+
  scale_fill_manual(values = my.cols)+
  theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                        axis.text.x  = element_text(size=20))+
  guides(fill = guide_legend(title="Site"))+ggtitle('Bandwidth (Hz)')+
  theme(legend.text=element_text(size=18))+
  theme(legend.title =element_text(size=20))+
  theme(title =element_text(size=20))+
  theme(axis.title.x =element_text(size=24, face="bold"))


site.intercepts <- extract(mfinal.stan, pars="site_rand_intercept")$site_rand_intercept
note.rate.intercepts <- data.frame(site.intercepts[, , 2]) # Intercepts for the 8th feature. 
names(note.rate.intercepts) <-levels(as.factor(combined.codas.all.sites$site))
str(note.rate.intercepts)

## Convert data into dataframe to pass to ggplot
CR <- cbind.data.frame(max.freq.intercepts$CR, rep("Crocker (CR)",length(max.freq.intercepts$CR)))
colnames(CR) <- c("samples","site")
DK <- cbind.data.frame(note.rate.intercepts$DK, rep("Deramakot (DK)",length(note.rate.intercepts$DK)))
colnames(DK) <- c("samples","site")
DV <- cbind.data.frame(note.rate.intercepts$DV, rep("Danum (DV)",length(note.rate.intercepts$DV)))
colnames(DV) <- c("samples","site")
IC <- cbind.data.frame(note.rate.intercepts$IC, rep("Imbak (IK)",length(note.rate.intercepts$IC)))
colnames(IC) <- c("samples","site")
KB <- cbind.data.frame(note.rate.intercepts$KB, rep("Kinabatangan (KB)",length(note.rate.intercepts$KB)))
colnames(KB) <- c("samples","site")
MB <- cbind.data.frame(note.rate.intercepts$MB, rep("Maliau (MB)",length(note.rate.intercepts$MB)))
colnames(MB) <- c("samples","site")
SAF <- cbind.data.frame(note.rate.intercepts$SA, rep("Kalabakan (KL)",length(note.rate.intercepts$SA)))
colnames(SAF) <- c("samples","site")
random.intercept.df.noterate <- rbind.data.frame(DK, DV, IC, KB,MB, SAF)

# Check structure of data frame
head(random.intercept.df.noterate)
dim(random.intercept.df.noterate)

## Create color palette for figure 
my.cols <- matlab::jet.colors(7)


## Create Figure 5 in ggplot
ggplot(random.intercept.df.noterate, aes(x=samples, fill=site))+ geom_density(alpha=.8, bw=.015)+
  xlab("Deviation from Base Line")+ylab("")+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                                  axis.text.x  = element_text(size=20))+
  #scale_fill_brewer(palette="OrRd")+
  scale_fill_manual(values = my.cols)+
  theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                        axis.text.x  = element_text(size=20))+
  guides(fill = guide_legend(title="Site"))+ggtitle('Note rate')+
  theme(legend.text=element_text(size=18))+
  theme(legend.title =element_text(size=20))+
  theme(title =element_text(size=20))+
  theme(axis.title.x =element_text(size=24, face="bold"))
