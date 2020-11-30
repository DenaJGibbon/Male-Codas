library(stringr)
library(tuneR)
library(seewave)
library(signal)

# Setting for MFCCs
min.freq = 500 
max.freq = 1500
n.windows = 9
num.cep = 12

files <- list.files('DuetTimingSelectionTables',recursive = T,
                    full.names = T,pattern = '.txt')

short.files <- list.files('DuetTimingSelectionTables',recursive = T,
                          full.names = F,pattern = '.txt')

short.files <- str_split_fixed(short.files,pattern = '/',n=2)[,2]

male.gps.data <- read.csv('male.coda.gps.csv')

all.mfcc.combined <- data.frame()
for(a in 1:length(files)){tryCatch({
  
  temp.file.path <- files[a]
  temp.table <- read.delim2(temp.file.path,stringsAsFactors = F)
  temp.table <- subset(temp.table,View=="Spectrogram 1")
  
  temp.table <- temp.table[order( as.numeric(temp.table$Begin.Time..s.)),]
  
  temp.table[,4:11] <- 
    temp.table[,4:11] %>% mutate_if(is.character,as.numeric)
  
  group <- str_split_fixed(short.files[a],pattern = '_',n=2)[,1]
  
  temp.wav <- str_split_fixed(short.files[a],pattern = '_',n=2)[,2]
  temp.wav.updated <- str_split_fixed(temp.wav,pattern = '[.]',n=2)[,1]
    
  temp.recording.info <- subset(male.gps.data,pair==group)

  wave.full.path <- paste(temp.recording.info$file.path,'/',group,'/',temp.wav.updated, '.WAV',sep='')
  wav.full <- tuneR::readWave(wave.full.path)
  
  list.sub <- which(temp.table$Call.type =='male.coda')
  
  if(length(list.sub)>0){
    
    list.codas <- split(list.sub, cumsum(c(1, diff(list.sub)) != 1))
    
    mfcc.vector.list <- list()
    for(c in 1:length(list.codas)){
      temp.coda.index <- list.codas[[c]]
      
      temp.table.coda <-  temp.table[min(temp.coda.index):max(temp.coda.index),]
      coda.start <- min(temp.table.coda$Begin.Time..s.) - 0.1
      coda.end <- max(temp.table.coda$End.Time..s.) + 0.1
      
      short.wav <- cutw(wav.full,from =coda.start , to=coda.end,output = 'Wave')
      
      # png(filename = paste('CodaSpectrograms/',group,'_',temp.wav.updated, c,'coda.png'), width=1000)
      # 
      # 
      # zcolors = colorRampPalette (c('dark blue','blue','cyan','light green','yellow',
      #                               'orange','red', 'brown'))
      # dynamicrange = 30
      # zrange = c(dynamicrange,0)
      # nlevels = abs (zrange[1] - zrange[2]) * 1.2
      # 
      # levels = pretty(zrange, nlevels)
      # zcolors = zcolors(length(levels) - 1)
      # 
      # Fs <-short.wav@samp.rate
      # step <- trunc(5*Fs/1000)             # one spectral slice every 5 ms
      # window <- trunc(20*Fs/1000)          # 40 ms data window
      # fftn <- 2^ceiling(log2(abs(window))) # next highest power of 2
      # spg <- specgram(short.wav@left, fftn, Fs, window, window-step)
      # S <- abs(spg$S[2:(fftn*3000/Fs),])   # magnitude in range 0<f<=4000 Hz.
      # S <- S/max(S)         # normalize magnitude so that max is 0 dB.
      # 
      # Sdb <- 20*log10(S) # Convert to dB
      # 
      # Sdb[which(Sdb < (-1 * dynamicrange))] = -1 * dynamicrange
      # 
      # image(t(Sdb),col=zcolors, axes = FALSE,useRaster = TRUE)
      # 
      # graphics.off()
      
      wav.dur <- duration(short.wav)
      win.time <- wav.dur/n.windows
      
      # Calculate MFCCs
      melfcc.output <- tuneR::melfcc(short.wav, minfreq = min.freq,
                                     hoptime = win.time, maxfreq = max.freq,
                                     numcep = num.cep, wintime = win.time)
      
      # Calculate delta cepstral coefficients
      deltas.output <- deltas(melfcc.output)
      
      # Ensure only same number of time windows are used for MFCC and delta coefficients Also append .wav duration
      mfcc.vector <- c(as.vector(t(melfcc.output[1:(n.windows - 1), 2:num.cep])), as.vector(t(deltas.output[1:(n.windows - 1), 2:num.cep])), wav.dur)
      
      mfcc.vector.list[[c]] <- mfcc.vector
    }
    
    
    mfcc.vector.df <- do.call(rbind.data.frame,mfcc.vector.list)
    colnames(mfcc.vector.df) <- seq(1,ncol(mfcc.vector.df),1)
    temp.mfcc.df <- cbind.data.frame(group, temp.wav.updated, mfcc.vector.df)
    all.mfcc.combined <- rbind.data.frame(all.mfcc.combined,temp.mfcc.df)
    write.csv(all.mfcc.combined,'all.mfcc.combined.csv',row.names = F)
  }
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

all.mfcc.combined.sub <- subset(all.mfcc.combined,X177 >=1)
all.mfcc.combined.sub$site <- substr(all.mfcc.combined.sub$group,start=1,stop=2)

all.mfcc.combined.sub <- na.omit(all.mfcc.combined.sub)

ml.model.svm <- e1071::svm(all.mfcc.combined.sub[,-c(1,2,180)], 
                           all.mfcc.combined.sub$group, kernel = "radial", 
                           cross = nrow(all.mfcc.combined.sub))

length(unique(all.mfcc.combined.sub$group))

ml.model.svm$tot.accuracy


# Assign to new object to modify


all.mfcc.combined.sub <- transform(all.mfcc.combined.sub,
                                      site=revalue(site,c("VJ"="SA",
                                                          "SC"="SA",
                                                          "SF"="SA",
                                                          "SM"="SA",
                                                          "CH"="SA")))

all.mfcc.combined.sub$site[which(all.mfcc.combined.sub$group=='CRIP')] <- 'SA'


all.mfcc.combined.sub$site <- as.factor(all.mfcc.combined.sub$site)

ml.model.svm.site <- e1071::svm(all.mfcc.combined.sub[,-c(1,2,180)], 
                           all.mfcc.combined.sub$site, kernel = "radial", 
                           cross = 25)


ml.model.svm.site$tot.accuracy


MaleCoda.umap <- 
  umap::umap(all.mfcc.combined.sub[,-c(1,2,180)],
             labels=as.factor(all.mfcc.combined.sub$group),
             controlscale=TRUE,scale=3)

recorder.id.df.umap <-
  cbind.data.frame(
    MaleCoda.umap$layout[,1:2],
    as.factor(all.mfcc.combined.sub$group)
  )


colnames(recorder.id.df.umap) <-
  c("Comp.1", "Comp.2", "site")


my_plot_umap <-
  ggplot(data = recorder.id.df.umap, aes(
    x = Comp.1,
    y = Comp.2,
    colour = site
  )) +
  geom_point(size = 3) +
  scale_color_manual(values = matlab::jet.colors(length(unique( as.factor(recorder.id.df.umap$site)))) ) +
  theme_bw() +xlab('UMAP: Dim 1')+ylab('UMAP: Dim 2')+ggtitle('Males')+ theme(legend.position = "none")

print(my_plot_umap)
