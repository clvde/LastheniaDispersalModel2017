fileloc <- "/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/May 28th 2017 - Burn In/"

startgen <- 1
totgens <- 345
popsizes <- startgen:totgens
maxlocs <- startgen:totgens
meandispa <- startgen:totgens
meandispb <- startgen:totgens
meanenv <- startgen:totgens
meanneut <- startgen:totgens
vardispa <- startgen:totgens
vardispb <- startgen:totgens
varenv <- startgen:totgens
varneut <- startgen:totgens
apply(gendata[,5:14], 1, sum)
sum(gendata[1,5:14])
for (i in startgen:totgens){
  filename <- paste("Burn_In_Gen",i,".csv",sep="")
  gendata <- read.csv(paste(fileloc,filename,sep=""))
  popsizes[i] <- nrow(gendata)
  maxlocs[i] <- max(gendata$location)
  meandispa[i] <- mean(apply(gendata[,5:14], 1, sum))
  meandispb[i] <- mean(apply(gendata[,15:24], 1, sum))
  meanenv[i] <- mean(apply(gendata[,25:34], 1, sum))
  meanneut[i] <- mean(apply(gendata[,35:44], 1, sum))
  vardispa[i] <- var(apply(gendata[,5:14], 1, sum))
  vardispb[i] <- var(apply(gendata[,15:24], 1, sum))
  varenv[i] <- var(apply(gendata[,25:34], 1, sum))
  varneut[i] <- var(apply(gendata[,35:44], 1, sum))
}
plot(maxlocs)
plot(popsizes)
plot(vardispb)
plot(vardispa)
plot(varenv)
plot(varneut)

dispa_table <- matrix(data=0, nrow = sum(popsizes), ncol = 3)
dispb_table <- matrix(data=0, nrow = sum(popsizes), ncol = 3)
env_table <- matrix(data=0, nrow = sum(popsizes), ncol = 3)
neut_table <- matrix(data=0, nrow = sum(popsizes), ncol = 3)

rowend = 0
rowstart = 0
for (i in 1:20){
  rowstart = rowend+1
  rowend = rowend + popsizes[i]
  print(c(rowstart,rowend))
  print(c(popsizes[i],rowend-rowstart))
}

for (i in 1:totgens){
  rowstart = rowend+1
  rowend = rowend + popsizes[i]
  generation <- matrix(data = i, nrow = popsizes[i], ncol = 1)
  dispa_table_onegen <- cbind(generation,locations, dispa_list)
  dispb_table_onegen <- cbind(generation,locations, dispb_list)
  env_table_onegen <- cbind(generation,locations, env_list)
  neut_table_onegen <- cbind(generation,locations, neut_list)
  dispa_table[rowstart:rowend] <- dispa_table_onegen
}

## ------- 

filename <- paste("Burn_In_Gen",990,".csv",sep="")
gendata <- read.csv(paste(fileloc,filename,sep=""))

locs <- gendata$location
cutted <- cut(locs, seq(-100,100),labels = FALSE)-100
binned_pre <- cbind(cutted, locs)
binned <- as.data.frame(binned_pre[order(cutted),])
head(binned)
binned$cutted
as.data.frame(binned)

DispParamData <- data.frame(bin = -100:100, mean_a = NA, mean_b = NA)
colnames(DispParamData)
head(DispParamData)
DispParamData$bin[100+101]

disp_a_list <- c()
disp_b_list <- c()
nrow(gendata)
for (i in 1:nrow(gendata))
  for (j in c(-100,100)){
    if (binned$cutted == j){
      cat(disp_a_list,sum(gendata[i,5:14])) 
    }
  }
