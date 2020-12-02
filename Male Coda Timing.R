# Differences in timing between great call and male coda?
library(stringr)
library(dplyr)
library(ggplot2)
library(lme4)
library(bbmle)
library(ggpubr)
library(plyr)
library(corrplot)
library(bayesplot)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())
library(stringr)

# Isolate file names
files <- list.files('DuetTimingSelectionTables',recursive = T,
                    full.names = T,pattern = '.txt')
short.files <- list.files('DuetTimingSelectionTables',recursive = T,
                          full.names = F,pattern = '.txt')

short.files <- str_split_fixed(short.files,pattern = '/',n=2)[,2]
  
# Omit based on multiple recording days
Files.ignore <- c(3,17,18,19,20,30,42,59,60,67,68,69,70,71,78,79,80)

files <- files[-Files.ignore]
short.files <- short.files[-Files.ignore]

# Read in selection tables
combined.coda.df <- data.frame()
for(a in 1:length(files)) { tryCatch({
  print(a)
  temp.table <- read.delim2(files[a],stringsAsFactors = F)
  temp.table <- temp.table[order( as.numeric(temp.table$Begin.Time..s.)),]
  
  group <- str_split_fixed(short.files[a],pattern = '_',n=2)[,1]
  group.label <- paste(group,a,sep='_')
  new.temp.table <-cbind.data.frame(temp.table,group.label)
  
  combined.coda.df <- rbind.data.frame(combined.coda.df,new.temp.table)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}




# solution
combined.coda.df[,4:11] <- 
  combined.coda.df[,4:11] %>% mutate_if(is.character,as.numeric)

pair.index <- unique(combined.coda.df$group.label)

pair.index <- pair.index[-c(35)]

coda.timing.df  <- data.frame()
for(b in 1:length(pair.index)){
  
  temp.df <- subset(combined.coda.df, group.label==pair.index[b])
  
  list.sub <- which(temp.df$Call.type =='great.call')
  
  if(length(list.sub)>0){
    
    list.codas <- split(list.sub, cumsum(c(1, diff(list.sub)) != 1))
    
    for(c in 1:length(list.codas)){
      
     great.call.end <-  temp.df[list.codas[[c]],]$End.Time..s.
     great.call.start <- temp.df[list.codas[[c]],]$Begin.Time..s.
     great.call.dur <-  great.call.end -great.call.start
     
     male.coda.temp <- temp.df[ (list.codas[[c]]+1),]
    
     
     if(c==length(list.codas)){
       male.coda.end <- temp.df[ nrow(temp.df),]
     } else{
       male.coda.end <- temp.df[ (list.codas[[c+1]])-1,]
     }
     
     if(nrow(male.coda.end) >1){
       
       male.coda.end <-  male.coda.end[which(male.coda.end$Call.type=='male.coda'),]
     }
     
     if( is.na(male.coda.temp$Call.type) != TRUE ){
     if(male.coda.temp$Call.type == 'male.coda' ){
      male.coda.start <-   male.coda.temp$Begin.Time..s.
    
     coda.timing <- male.coda.start-great.call.end
     coda.duration <- male.coda.end$End.Time..s. -male.coda.temp$Begin.Time..s.
     
     individual <- str_split_fixed( temp.df$group.label[1],pattern = '_',n=2)[,1]
     call.id <- unique(paste(temp.df$group.label,c,sep='.'))  
     
     coda.timing.temp <- cbind.data.frame(coda.timing,individual,call.id,coda.duration)
     coda.timing.df <- rbind.data.frame(coda.timing.df,coda.timing.temp)
      }
     }
    }
  }
}

coda.timing.df <- droplevels(subset(coda.timing.df,coda.timing <5))
coda.timing.df$site <- substr(coda.timing.df$individual,start=1,stop=2)
coda.timing.df  <- transform(coda.timing.df ,
                                      site=revalue(site,c("VJ"="SA",
                                                          "SC"="SA",
                                                          "SF"="SA",
                                                          "SM"="SA",
                                                          "CH"="SA")))

coda.timing.df$site[which(coda.timing.df $individual=='CRIP')] <- 'SA'


ggpubr::ggdensity(data=coda.timing.df,x="coda.timing",fill = "site")

coda.timing.df$individual <- as.factor(as.integer(coda.timing.df$individual))
ggpubr::ggboxplot(data=coda.timing.df,x="individual" , y="coda.timing",fill = "site",
                  palette = matlab::jet.colors(length(unique(coda.timing.df$site))))+
  xlab('Male')+ ylab('Coda timing relative to female (s)')+rotate_x_text(angle = 90)

ggpubr::ggerrorplot(data=coda.timing.df,x="individual" , y="coda.duration",color = "site")

hist(coda.timing.df$coda.duration)
hist(coda.timing.df$coda.timing)
table(coda.timing.df$individual)
length(unique(coda.timing.df$individual))
nrow(coda.timing.df)

# MANOVA for duet timing
## Isolate relevant features from data set
d.manova.timing <- coda.timing.df[,c("coda.timing","coda.duration")]


## Check the structure of the data
cor(d.manova.timing)

## Log transform data
d.manova.timing$coda.duration <- log(d.manova.timing$coda.duration)

pairs(d.manova.timing)

## Check the structure of the data
str(d.manova.timing)

### Set-up data to pass to Stan. 
# Integer-coded vector of group IDs
group.int <- as.numeric(coda.timing.df$individual)

# Check structure 
table(group.int)

# Integer-coded vector of site IDs
site.int <- as.numeric(as.factor(coda.timing.df$site))

# Check structure 
table(site.int)

# Center data matrix at feature means
col.means <- apply(d.manova.timing, MARGIN=2, FUN="mean")
y.centered <- sweep(d.manova.timing, MARGIN=2, STATS=col.means)

# Create a data list to pass to Stan
data_list <- list(
  K = dim(d.manova.timing)[2],
  J= length(unique(coda.timing.df $individual)),
  M= length(unique(coda.timing.df $site)),
  N= dim(d.manova.timing)[1],
  y= as.matrix(y.centered), 
  group= group.int,
  site= site.int
)


# Code to run the STAN model
# NOTE: the .stan file must be linked in the code below
coda.timing.stan = stan(file="MANOVA.male codas.stan", 
                   model_name = "mfinal", 
                   data=data_list, iter=3000, warmup=1500, chains=2, 
                   cores=2, 
                   control = list(stepsize = 0.5, adapt_delta = 0.99, max_treedepth = 20))

# Optional code to save the output
#save(mfinal.stan, file = "/Volumes/DJC HardDrive/stan.model.output.oct.2017.rda")

## Check model output
# Create traceplots to check for mixing for site level variance
draws <- as.array(coda.timing.stan, pars="ICC_site")
mcmc_trace(draws)
stan_dens(coda.timing.stan, pars=c("ICC_site"))
round(summary(coda.timing.stan, pars=c("ICC_site"))$summary, 3)

# Create traceplots to check for mixing for group level variance
draws <- as.array(coda.timing.stan, pars="ICC_group")
mcmc_trace(draws)
stan_dens(coda.timing.stan, pars=c("ICC_group"))
round(summary(coda.timing.stan, pars=c("ICC_group"))$summary, 3)

## Extract site-specific random intercepts for trill rate
site.intercepts <- extract(coda.timing.stan, pars="site_rand_intercept")$site_rand_intercept
trill.intercepts <- data.frame(site.intercepts[, , 1]) # Intercepts for the 8th feature. 
names(trill.intercepts) <-levels(as.factor(coda.timing.df $site))
str(trill.intercepts)

## Convert data into dataframe to pass to ggplot
DK <- cbind.data.frame(trill.intercepts$DK, rep("Deramakot (DK)",length(trill.intercepts$DK)))
colnames(DK) <- c("samples","site")
DV <- cbind.data.frame(trill.intercepts$DV, rep("Danum (DV)",length(trill.intercepts$DV)))
colnames(DV) <- c("samples","site")
IC <- cbind.data.frame(trill.intercepts$IC, rep("Imbak (IK)",length(trill.intercepts$IC)))
colnames(IC) <- c("samples","site")
MB <- cbind.data.frame(trill.intercepts$MB, rep("Maliau (MB)",length(trill.intercepts$MB)))
 colnames(MB) <- c("samples","site")
SAF <- cbind.data.frame(trill.intercepts$SA, rep("Kalabakan (KL)",length(trill.intercepts$SA)))
colnames(SAF) <- c("samples","site")
random.intercept.df <- rbind.data.frame(DK, DV, IC, SAF)

# Check structure of data frame
head(random.intercept.df)
dim(random.intercept.df)

## Create color palette for figure 
my.cols <- viridis::viridis(n=7, option="D")


## Create Figure 5 in ggplot
ggplot(random.intercept.df, aes(x=samples, fill=site))+ geom_density(alpha=.8, bw=.015)+
  xlab("Deviation from Base Line")+ylab("")+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                                                  axis.text.x  = element_text(size=20))+
  #scale_fill_brewer(palette="OrRd")+
  scale_fill_manual(values = my.cols)+
  theme_classic()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                        axis.text.x  = element_text(size=20))+
  guides(fill = guide_legend(title="Site"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.title =element_text(size=20))+
  theme(axis.title.x =element_text(size=24, face="bold"))

    