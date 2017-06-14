library(dplyr)
bi = 1

# Goes from Trait files to summary statistics for head trait in each generation 

for (g in seq(from=5100, to=10100, by=100)){
  # Specify file name
  print(g)
  traitfileloc <- paste("/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/No Change Extension ",bi," Trait Data/",sep="")
  traitfilename <- paste("Trait_Data_NoChange_Extension_",bi,"_Gen_",g,".csv",sep="")

  # Read in file
  traitdata <- read.csv(paste(traitfileloc, traitfilename, sep=""))

  # Grab location data, bin it, and stick bin column into data frame
  raw_locs <- traitdata$Locations
  BinnedLocations <- cut(raw_locs, seq(-99.5,100.5),labels = FALSE)-100
  traitdata_mod <- cbind(traitdata[,2:4],BinnedLocations,traitdata[,5:8])
  traitdata_mod <- traitdata_mod[order(traitdata_mod[,4]),]

  # create a list of the different bins (midpoints of bins that are 1 unit wide. So bin 0 goes from -0.5 to 0.5)
  uniq <- unique(traitdata_mod$BinnedLocations)

  # Initialize data frames to store lists of different phenotype data summarized within locations. One data frame for each phenotype

  Dispa_analyzed <- matrix(data = NA, nrow = (length(uniq)-1), ncol = 6)
  Dispb_analyzed <- matrix(data = NA, nrow = (length(uniq)-1), ncol = 6)
  Env_analyzed <- matrix(data = NA, nrow = (length(uniq)-1), ncol = 6)
  Neut_analyzed <- matrix(data = NA, nrow = (length(uniq)-1), ncol = 6)

  for (i in 1:length(uniq)){
    bin <- uniq[i]
    filty <- filter(traitdata_mod,traitdata_mod$BinnedLocations == bin)
    nobsv <- nrow(filty)
  
    meanvec <- as.vector(apply(filty[,5:8],2,mean))
    maxvec <- as.vector(apply(filty[,5:8],2,max))
    minvec <- as.vector(apply(filty[,5:8],2,min))
    varvec <- as.vector(apply(filty[,5:8],2,var))
  
    Dispa_analyzed[i,] <-  c(bin, meanvec[1],maxvec[1],minvec[1],varvec[1],nobsv)
    Dispb_analyzed[i,] <- c(bin, meanvec[2], maxvec[2], minvec[2], varvec[2], nobsv)
    Env_analyzed[i,] <- c(bin, meanvec[3], maxvec[3], minvec[3], varvec[3], nobsv)
    Neut_analyzed[i,] <- c(bin, meanvec[4], maxvec[4], minvec[4], varvec[4], nobsv)
  }

  cnms <- c("LocationBin", "Mean", "Maximum", "Minimum", "Variance", "NumObservations")

  Dispa_analyzed <- as.data.frame(Dispa_analyzed)
  Dispb_analyzed <- as.data.frame(Dispb_analyzed)
  Env_analyzed <- as.data.frame(Env_analyzed)
  Neut_analyzed <- as.data.frame(Neut_analyzed)
  
  colnames(Dispa_analyzed) <- cnms
  colnames(Dispb_analyzed) <- cnms
  colnames(Env_analyzed) <- cnms
  colnames(Neut_analyzed) <- cnms

  write.csv(Dispa_analyzed,paste(traitfileloc, "Dispa_analyzed_NoChange_Extension_",bi,"_Gen_",g,".csv", sep=""))
  write.csv(Dispb_analyzed,paste(traitfileloc, "Dispb_analyzed_NoChange_Extension_",bi,"_Gen_",g,".csv", sep=""))
  write.csv(Env_analyzed,paste(traitfileloc, "Env_analyzed_NoChange_Extension_",bi,"_Gen_",g,".csv", sep=""))
  write.csv(Neut_analyzed,paste(traitfileloc, "Neut_analyzed_NoChange_Extension_",bi,"_Gen_",g,".csv", sep=""))
  
}

# Visualizing Spread
ana2600 <- read.csv(paste(traitfileloc, "Dispa_analyzed_NoChange_Extension_",bi,"_Gen_",2600,".csv", sep=""))
gendata <- read.csv("/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/June 7th 2017 Burn In 2/Burn_In_Gen1100.csv")


## Trying to make data frames that will work well for making heatmaps!!

# define the largest generation number that has been completed. When analyzing multiple 
maxgen <- 10100

locseq <- seq(from=-400, to=400)
timeseq <- seq(from=100, to=10100)

# rows are time points (i.e. every 100 generations or something), cols are spatial locations.
heatmap_means <- matrix(nrow = 6, ncol = length(locseq))
heatmap_maxs <- matrix(nrow = 6, ncol = length(locseq))
heatmap_mins <- matrix(nrow = 6, ncol = length(locseq))
heatmap_vars <- matrix(nrow = 6, ncol = length(locseq))


traitfileloc <- paste("/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/No Change Extension ",bi," Trait Data/",sep="")
cf <- read.csv(paste(traitfileloc, "Dispa_analyzed_NoChange_Extension_",1,"_Gen_",5100,".csv", sep=""))
currentfile <- cf[,-1]

is.data.frame(currentfile)

for (i in seq(from=5100, to=10100, by=100)){
  currentfile <- read.csv(paste(traitfileloc, "Dispa_analyzed_NoChange_Extension_",bi,"_Gen_",i,".csv", sep=""))
  locate <- currentfile$LocationBin
  means <- currentfile$Mean
  # for (j in locate){
  #   cn <- paste(j,sep="")
  #   print(heatdf_mean$cn)
  # }
}

# Plotting!
par(mfrow=c(2,2))
plot(x=Dispa_analyzed$LocationBin, y=Dispa_analyzed$Mean,ylim=c(0.2,0.3))
plot(x=Dispb_analyzed$LocationBin, y=Dispb_analyzed$Mean,ylim=c(0.95,1.05))
plot(x=Env_analyzed$LocationBin, y=Env_analyzed$Mean,ylim=c(-1.1,-0.9))
plot(x=Neut_analyzed$LocationBin, y=Neut_analyzed$Mean,ylim=c(-0.1,0.1))

par(mfrow=c(2,2))
plot(x=Dispa_analyzed$LocationBin, y=Dispa_analyzed$Maximum,ylim=c(0.2,0.3))
plot(x=Dispb_analyzed$LocationBin, y=Dispb_analyzed$Maximum,ylim=c(0.95,1.05))
plot(x=Env_analyzed$LocationBin, y=Env_analyzed$Maximum,ylim=c(-1.1,-0.9))
plot(x=Neut_analyzed$LocationBin, y=Neut_analyzed$Maximum,ylim=c(-0.1,0.1))

par(mfrow=c(2,2))
plot(x=Dispa_analyzed$LocationBin, y=Dispa_analyzed$Minimum,ylim=c(0.2,0.3))
plot(x=Dispb_analyzed$LocationBin, y=Dispb_analyzed$Minimum,ylim=c(0.95,1.05))
plot(x=Env_analyzed$LocationBin, y=Env_analyzed$Minimum,ylim=c(-1.1,-0.9))
plot(x=Neut_analyzed$LocationBin, y=Neut_analyzed$Minimum,ylim=c(-0.1,0.1))

par(mfrow=c(2,2))
plot(x=Dispa_analyzed$LocationBin, y=Dispa_analyzed$Variance,ylim=c(0,0.01))
plot(x=Dispb_analyzed$LocationBin, y=Dispb_analyzed$Variance,ylim=c(0,0.01))
plot(x=Env_analyzed$LocationBin, y=Env_analyzed$Variance,ylim=c(0,0.01))
plot(x=Neut_analyzed$LocationBin, y=Neut_analyzed$Variance,ylim=c(0,0.01))