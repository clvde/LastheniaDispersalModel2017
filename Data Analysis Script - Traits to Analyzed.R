library(dplyr)
bi = 3

# Goes from Trait files to summary statistics for head trait in each generation 

for (g in seq(from=5100, to=30000, by=100)){
  print(g)
  
  # Specify file name
  traitfileloc <- paste("/Users/Courtney/Documents/Melbourne & Flaxman Labs/Simulation Results/Env Change ",bi," (v=0.01) Trait Data/",sep="")
  traitfilename <- paste("Trait_Data_Env_Change_(v=0.01)_",bi,"_Gen_",g,".csv",sep="")

  # Read in file
  traitdata <- read.csv(paste(traitfileloc, traitfilename, sep=""))

  # Grab location data, bin it, and stick bin column into data frame
  raw_locs <- traitdata$Locations
  BinnedLocations <- cut(raw_locs, seq(-700.5,700.5), labels = FALSE)-701
  traitdata_mod <- cbind(traitdata[,2:4],BinnedLocations,traitdata[,5:8])
  traitdata_mod <- traitdata_mod[order(traitdata_mod[,4]),]
  
  # create a list of the different bins (midpoints of bins that are 1 unit wide. So bin 0 goes from -0.5 to 0.5)
  uniq <- unique(traitdata_mod$BinnedLocations)
  
  # Initialize data frames to store lists of different phenotype data summarized within locations. One data frame for each phenotype
  
  Dispa_analyzed <- matrix(data = NA, nrow = length(uniq), ncol = 6)
  Dispb_analyzed <- matrix(data = NA, nrow = length(uniq), ncol = 6)
  Env_analyzed <- matrix(data = NA, nrow = length(uniq), ncol = 6)
  Neut_analyzed <- matrix(data = NA, nrow = length(uniq), ncol = 6)
  
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

  write.csv(Dispa_analyzed,paste(traitfileloc, "Dispa_analyzed_Env_Change_(v=0.01)_",bi,"_Gen_",g,".csv", sep=""))
  write.csv(Dispb_analyzed,paste(traitfileloc, "Dispb_analyzed_Env_Change_(v=0.01)_",bi,"_Gen_",g,".csv", sep=""))
  write.csv(Env_analyzed,paste(traitfileloc, "Env_analyzed_Env_Change_(v=0.01)_",bi,"_Gen_",g,".csv", sep=""))
  write.csv(Neut_analyzed,paste(traitfileloc, "Neut_analyzed_Env_Change_(v=0.01)_",bi,"_Gen_",g,".csv", sep=""))
  
}