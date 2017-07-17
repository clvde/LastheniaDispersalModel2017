## Plotting Genetic Variation

## -------------
## Burn In 
bi <- 3
fileloc <- paste("/Volumes/Seagate Backup Plus Drive/Dispersal Evolution Model Data/Env Change ",bi," (v=0.01)/",sep="")
maxgen <- 3950
mingen <- 10
totlen <- length(seq(from=mingen, to=maxgen, by=10))

dispa_varlist <- vector(length=totlen)
dispb_varlist <- vector(length=totlen)
env_varlist <- vector(length=totlen)
neut_varlist <- vector(length=totlen)

for (i in seq(from=mingen, to=5000, by=10)){
  filename <- paste("Burn_In",bi,"_Gen",i,".csv",sep="")
  gendata <- read.csv(paste(fileloc,filename,sep=""))
  LocCorrData <- cbind(gendata$generation+1, gendata$individual_ID, gendata$location, apply(gendata[,5:14], 1, sum), apply(gendata[,15:24], 1, sum), apply(gendata[,25:44], 1, sum), apply(gendata[,45:64], 1, sum))
  LocCorrData <- as.data.frame(LocCorrData)
  colnames(LocCorrData) <- c("Generation", "Individual_ID", "Locations", "Disp_a_trait", "Disp_b_trait", "Env_trait", "Neut_trait")
  dispa_varlist[(i/10)] = var(LocCorrData$Disp_a_trait)
  dispb_varlist[(i/10)] = var(LocCorrData$Disp_b_trait)
  env_varlist[(i/10)] = var(LocCorrData$Env_trait)
  neut_varlist[(i/10)] = var(LocCorrData$Neut_trait)
  write.csv(LocCorrData, paste("/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/Burn In ",bi," Trait Data/Trait_Data_Burn_In_",bi,"_Gen_",i,".csv",sep=""))
} 

## -------------
## Extension

bi <- 3
maxgen <- 25000
mingen <- 5010
totlen <- length(seq(from=mingen, to=maxgen, by=10))

fileloc <- paste("/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/No Change Extension ",bi,"/",sep="")
 
for (i in seq(from=mingen, to=maxgen, by=10)){
  print(i)
  filename <- paste("NoChange_Extension_Gen",i,".csv",sep="")
  gendata <- read.csv(paste(fileloc,filename,sep=""))
  LocCorrData <- cbind(gendata$generation+1, gendata$individual_ID, gendata$location, apply(gendata[,5:14], 1, sum), apply(gendata[,15:24], 1, sum), apply(gendata[,25:44], 1, sum), apply(gendata[,45:64], 1, sum))
  LocCorrData <- as.data.frame(LocCorrData)
  colnames(LocCorrData) <- c("Generation", "Individual_ID", "Locations", "Disp_a_trait", "Disp_b_trait", "Env_trait", "Neut_trait")
  write.csv(LocCorrData, paste("/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/No Change Extension ",bi," Trait Data/Trait_Data_NoChange_Extension_",bi,"_Gen_",i,".csv",sep=""))
}

## -------------
## Env Change

bi <- 3
fileloc <- paste("/Volumes/Seagate Backup Plus Drive/Dispersal Evolution Model Data/Env Change ",bi," (v=0.01)/",sep="")
maxgen <- 25000
mingen <- 10
totlen <- length(seq(from=mingen, to=maxgen, by=10))

dispa_varlist <- vector(length=totlen)
dispb_varlist <- vector(length=totlen)
env_varlist <- vector(length=totlen)
neut_varlist <- vector(length=totlen)

for (i in seq(from=mingen, to=maxgen, by=10)){
  print(i)
  filename <- paste("Env_Change_F_BI",bi,"_Gen",i,".csv",sep="")
  gendata <- read.csv(paste(fileloc,filename,sep=""))
  LocCorrData <- cbind(gendata$generation+5001, gendata$individual_ID, gendata$location, apply(gendata[,5:14], 1, sum), apply(gendata[,15:24], 1, sum), apply(gendata[,25:44], 1, sum), apply(gendata[,45:64], 1, sum))
  LocCorrData <- as.data.frame(LocCorrData)
  colnames(LocCorrData) <- c("Generation", "Individual_ID", "Locations", "Disp_a_trait", "Disp_b_trait", "Env_trait", "Neut_trait")
  #dispa_varlist[(i/10)] = var(LocCorrData$Disp_a_trait)
  #dispb_varlist[(i/10)] = var(LocCorrData$Disp_b_trait)
  #env_varlist[(i/10)] = var(LocCorrData$Env_trait)
  #neut_varlist[(i/10)] = var(LocCorrData$Neut_trait)
  write.csv(LocCorrData, paste("/Users/Courtney/Documents/Melbourne & Flaxman Labs/Simulation Results/Env Change ",bi," (v=0.01) Trait Data/Trait_Data_Env_Change_(v=0.01)_",bi,"_Gen_",(i+5000),".csv",sep=""))
}

par(mfrow=c(2,2))
plot(dispa_varlist)
plot(dispb_varlist)
plot(env_varlist)
plot(neut_varlist)

## -----------------------------------------------
## Change file names on no change extension 1 data

fileloc <- "/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/No Change Extension 3/"

for (i in seq(from=14590, to=16000, by=10)){
  filename <- paste("No_Change_BI3_Gen",i,".csv",sep="")
  gendata <- read.csv(paste(fileloc,filename,sep=""))
  
  c_p2 <- gendata[,-1]
  c_p2[,1] <- c_p2[,1]+1
  current_population <- c_p2
  print(paste("NoChange_Extension3_Gen",i,".csv",sep=""))
  write.csv(current_population, paste(fileloc,"NoChange_Extension3_Gen",i,".csv", sep=""))
}

