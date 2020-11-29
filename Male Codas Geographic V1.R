library(stringr)
library(dplyr)
library(ggplot2)
library(lme4)
library(bbmle)
library(plyr)

# V1 focus only on codas that were annotated in the same way

files <- list.files('DuetTimingSelectionTables',recursive = T,
           full.names = T,pattern = '.txt')
short.files <- list.files('DuetTimingSelectionTables',recursive = T,
                          full.names = F,pattern = '.txt')

short.files <- str_split_fixed(short.files,pattern = '/',n=2)[,2]

# Omit based on multiple recording days
Files.ignore <- c(3,17,18,19,20,30,42,59,60,67,68,69,70,71,78,79,80)

files <- files[-Files.ignore]
short.files <- short.files[-Files.ignore]

# Code to read in selection tables
combined.coda.df <- data.frame()
for(a in 1:length(files)) { tryCatch({
  print(a)
  temp.table <- read.delim2(files[a],stringsAsFactors = F)
  temp.table <- subset(temp.table,View=="Spectrogram 1")
  
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

malecodadf <- data.frame()
for(b in 1:length(pair.index)){

temp.df <- subset(combined.coda.df, group.label==pair.index[b])

list.sub <- which(temp.df$Call.type =='male.coda')

if(length(list.sub)>0){

list.codas <- split(list.sub, cumsum(c(1, diff(list.sub)) != 1))

for(c in 1:length(list.codas)){
  temp.coda.index <- list.codas[[c]]
  
  temp.table <-  temp.df[min(temp.coda.index):max(temp.coda.index),]
individual <- str_split_fixed( temp.table$group.label[1],pattern = '_',n=2)[,1]
note.dur <- mean(temp.table$End.Time..s. - temp.table$Begin.Time..s.)  
call.dur  <- max(temp.table$End.Time..s.) - min(temp.table$Begin.Time..s.)
call.id <- unique(paste(temp.table$group.label,c,sep='.'))

mean5 <- mean(temp.table$Freq.5...Hz.)
mean95 <-  mean(temp.table$Freq.95...Hz.)
max5 <- max(temp.table$Freq.5...Hz.)
max95 <-  max(temp.table$Freq.95...Hz.)
min5 <- min(temp.table$Freq.5...Hz.)
min95 <-  min(temp.table$Freq.95...Hz.)
minbw <-  min(temp.table$BW.90...Hz.)
maxbw <-  max(temp.table$BW.90...Hz.)
meanbw <- mean(temp.table$BW.90...Hz.)
mindurnote <-  min(temp.table$Dur.90...s.)
maxdurnote <- max(temp.table$Dur.90...s.)
nnotes <- nrow(temp.table)
noterate <- nrow(temp.table)/call.dur
note1dur <- temp.table[1,]$Dur.90...s.
note1minfreq <- temp.table[1,]$Freq.5...Hz.
note1maxfreq <- temp.table[1,]$Freq.95...Hz.
range.bw <- maxbw-minbw
rest.dur <- call.dur - sum(temp.table$Dur.90...s.)
lastnotedur <- temp.table[nrow(temp.table),]$Dur.90...s.
lastnoteminfreq <- temp.table[nrow(temp.table),]$Freq.5...Hz.
lastnotemaxfreq <- temp.table[nrow(temp.table),]$Freq.95...Hz.

temp.coda.df <- cbind.data.frame(individual,call.id,call.dur,nnotes,
                                  min5,min95, minbw,maxbw,mean5,mean95,max5,max95,meanbw,mindurnote,
                                 maxdurnote,noterate,note1dur,note1minfreq,note1maxfreq,range.bw,rest.dur,
                                 lastnotedur,lastnoteminfreq,lastnotemaxfreq )

malecodadf <- rbind.data.frame(malecodadf,temp.coda.df)
}
}
}

malecodadf$site <- substr(malecodadf$individual,start=1,stop=2)

# Assign to new object to modify
combined.codas.all.sites <- malecodadf

combined.codas.all.sites <- transform(combined.codas.all.sites,
                            site=revalue(site,c("VJ"="SA",
                                                "SC"="SA",
                                                "SF"="SA",
                                                "SM"="SA",
                                                "CH"="SA")))

combined.codas.all.sites$site[which(combined.codas.all.sites$individual=='CRIP')] <- 'SA'

combined.codas.all.sites$Site <- as.factor(combined.codas.all.sites$site)


# See increase in complexity over course of duet?
pair.id <- unique(combined.codas.all.sites$individual)

complexity.df <- data.frame()
for(d in 1:length(pair.id)){
  tempdf <- subset(combined.codas.all.sites,individual==pair.id[d])
  tempdf$seq <- seq(1,nrow(tempdf),1)
  complexity.df <- rbind.data.frame(complexity.df,tempdf)
}

# Complexity as number of notes
complexitymodel <- glmer(nnotes ~ seq + (1|site/individual), family='poisson',  data=complexity.df)
complexitymodel.null <- glmer(nnotes ~ 1 +  (1|site/individual), family='poisson',  data=complexity.df)

AICctab(complexitymodel,complexitymodel.null)

coefplot::coefplot(complexitymodel,intercept=F)

sjPlot::plot_model(complexitymodel,type='eff')

MuMIn::r.squaredGLMM(complexitymodel)

# Note rate?
complexitymodel.noterate <- lmer(noterate ~ seq + (1|site/individual),   data=complexity.df)
complexitymodel.noterate.null <- lmer(noterate ~ 1 +  (1|site/individual),   data=complexity.df)

AICctab(complexitymodel.noterate,complexitymodel.noterate.null)

coefplot::coefplot(complexitymodel.noterate,intercept=F)

sjPlot::plot_model(complexitymodel.noterate.null,type='re')

# Coda duration
complexitymodel.call.dur <- lmer(call.dur ~ seq + (1|site/individual),   data=complexity.df)
complexitymodel.call.dur.null <- lmer(call.dur ~ 1 +  (1|site/individual),   data=complexity.df)

AICctab(complexitymodel.call.dur,complexitymodel.call.dur.null)

coefplot::coefplot(complexitymodel.call.dur,intercept=F)

sjPlot::plot_model(complexitymodel.call.dur,type='eff')

ggpubr::ggdensity(data=complexity.df,x='call.dur',fill='site')

coefplot::multiplot(complexitymodel,complexitymodel.call.dur,
                    intercept = F)+theme_bw()

# Male coda types
# cluster.df <- apcluster::apcluster(
#   negDistMat(r = 2),
#   combined.codas.all.sites[,c(3:19)],
#   maxits = 10000,
#   convits = 10000,
#   nonoise = T,
#   q=0.1
# )
# 
# length(cluster.df@exemplars)
# 
# complexity.df$cluster <- as.factor(unlist(cluster.df@idx))
# 
# cluster.model <- lmer(cluster ~ seq + (1|site/individual),   data=complexity.df)
# 
# b1 <- brm (cluster ~ seq + + (1|site/individual),
#            data=complexity.df, family="categorical",
#            prior=c(set_prior ("normal (0, 8)")))


## UMAP

male.umap <- 
  umap::umap(combined.codas.all.sites[,c(3:13,16)],labels=as.factor(combined.codas.all.sites$Site),
             controlscale=TRUE,scale=3)

plot.for.males <-
  cbind.data.frame(male.umap$layout[,1:2],
                   combined.codas.all.sites$Site)
colnames(plot.for.males) <-
  c("Dim.1", "Dim.2", "Site")

my_plot_males <-
  ggplot(data = plot.for.males, aes(
    x = Dim.1,
    y = Dim.2,
    colour = Site
  )) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(plot.for.males$Site)))) +
  theme_bw() + ggtitle('Male codas') + xlab('UMAP: Dim 1')+ylab('UMAP: Dim 2')

my_plot_males


ml.model.svm.site <- e1071::svm(combined.codas.all.sites[,c(3:24)], 
                                combined.codas.all.sites$Site, kernel = "radial", 
                                cross = 25)


ml.model.svm.site$tot.accuracy

# Individual
combined.codas.all.sites$individual <- as.factor(combined.codas.all.sites$individual)

male.individual.umap <- 
  umap::umap(combined.codas.all.sites[,c(3:13,16)],labels=as.factor(combined.codas.all.sites$individual),
             controlscale=TRUE,scale=3)

plot.for.male.individuals <-
  cbind.data.frame(male.individual.umap$layout[,1:2],
                   combined.codas.all.sites$individual)
colnames(plot.for.male.individuals) <-
  c("Dim.1", "Dim.2", "individual")

my_plot_male.individuals <-
  ggplot(data = plot.for.male.individuals, aes(
    x = Dim.1,
    y = Dim.2,
    colour = individual
  )) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors (length(unique(plot.for.male.individuals$individual)))) +
  theme_bw() + ggtitle('Male codas') + xlab('UMAP: Dim 1')+ylab('UMAP: Dim 2')+ theme(legend.position = "none")

my_plot_male.individuals


ml.model.svm <- e1071::svm(combined.codas.all.sites[,c(3:24)], 
                           combined.codas.all.sites$individual, kernel = "radial", 
                           cross = 25)


ml.model.svm$tot.accuracy

