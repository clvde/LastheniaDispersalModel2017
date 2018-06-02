library(ggplot2)

#------------------
# CREATE DATA FRAME
#------------------

## Another way to organize the data for production of a heatmap

bi <- 3
maxgen_bi <- 5000
mingen_bi <- 100
maxgen <- 30000
mingen <- 5100
trait <- "Env"

traitfileloc <- paste("/Users/Courtney/Documents/Melbourne & Flaxman Labs/Simulation Results/Burn In ",bi," Trait Data/",sep="")
currentfile <- read.csv(paste(traitfileloc, trait, "_analyzed_Burn_In_",bi,"_Gen_",mingen_bi,".csv", sep=""))
currentfile <- currentfile[,-1]
currentfile <- cbind(mingen,currentfile)
bigdf <- currentfile
colnames(bigdf)[1] <- "Generation"

for (i in seq(from=mingen_bi+100, to=maxgen_bi, by=100)){
  currentfile <- read.csv(paste(traitfileloc,trait,"_analyzed_Burn_In_",bi,"_Gen_",i,".csv", sep=""))
  currentfile <- currentfile[,-1]
  currentfile <- cbind(i,currentfile)
  colnames(currentfile)[1] <- "Generation"
  
  bigdf <- rbind(bigdf,currentfile)
}
# 
traitfileloc <- paste("/Users/Courtney/Documents/Melbourne & Flaxman Labs/Simulation Results/Env Change ",bi," (v=0.01) Trait Data/",sep="")
 
for (i in seq(from=mingen, to=maxgen, by=100)){
  currentfile <- read.csv(paste(traitfileloc,trait,"_analyzed_Env_Change_(v=0.01)_",bi,"_Gen_",i,".csv", sep=""))
  currentfile <- currentfile[,-1]
  currentfile <- cbind(i,currentfile)
  colnames(currentfile)[1] <- "Generation"
 
  bigdf <- rbind(bigdf,currentfile)
}
write.csv(bigdf,paste(traitfileloc,trait,"_BIandEnvChange(v=0.01)_",bi,"_for_heatmap.csv",sep=""))

bigdfsub <- bigdf[-dim(bigdf)[1],]
bigdfsub[116956:116958,]
#---------------------------
# PLOT: heatmap
# - here, we use geom_tile()
#---------------------------

#Dispa
min(bigdfsub$Max)
max(bigdfsub$Max)
b <- c(-1.06,-1,-0.94)
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + scale_fill_gradientn(limits = c(-1.06,-0.94), colours=c("#0072ff", "white", "#ff0000"), breaks=b, labels=format(b)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "black"), axis.line = element_line(colour = "black")) + geom_raster(aes(fill = Mean))
                                                                                                                                                                                                                                                                                                                                                                                                

#Dispb
max(bigdfsub$Mean)
min(bigdfsub$Mean)
b <- c(0.94,1,1.06)
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + scale_fill_gradientn(limits = c(0.94,1.06), colours=c("#0090ff", "white", "darkred"), breaks=b, labels=format(b)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "black"), axis.line = element_line(colour = "black")) + geom_raster(aes(fill = Mean), interpolate = TRUE)

#Env
max(bigdfsub$Mean)
min(bigdfsub$Mean)
b <- c(-1.07,-1,-0.93)
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + scale_fill_gradientn(limits = c(-1.07,-0.93), colours=c("#0090ff", "white", "darkred"), breaks=b, labels=format(b)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = '#919191', colour = 'white'), axis.line = element_line(colour = "black")) + geom_tile(aes(fill = Mean))

#Neut
max(bigdfsub$Mean)
min(bigdfsub$Mean)
b <- c(-0.07,0,0.07)
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + scale_fill_gradientn(limits = c(-0.07,0.07), colours=c("#0090ff", "white", "darkred"), breaks=b, labels=format(b)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = '#919191', colour = 'white'), axis.line = element_line(colour = "black")) + geom_tile(aes(fill = Mean))

## -----
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + geom_tile(aes(fill = Mean))


sort(bigdfsub$Variance)[77000:78049]
length(bigdfsub$Variance)
b <- c(0,0.0006)
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + scale_fill_gradientn(limits = c(0,0.0006), colours=c("black", "white"), breaks=b, labels=format(b)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_tile(aes(fill = Variance))




ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + geom_tile(aes(fill = Maximum)) 
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + geom_tile(aes(fill = Minimum)) 
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + geom_tile(aes(fill = Variance))
ggplot(data = bigdfsub, aes(x = LocationBin, y = Generation)) + geom_tile(aes(fill = NumObservations)) 
