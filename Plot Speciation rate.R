## This is an experimental code that would need some work to be more user-friendly

if ( ! ("picante" %in% installed.packages())) {install.packages("picante", dependencies=T)}
if ( ! ("pspline" %in% installed.packages())) {install.packages("pspline", dependencies=T)}

library("picante")
library("pspline")

pdf("Rates as a function of paleo-envrionmental variables.pdf")

par(mfrow=c(2,3), mar=c(4,4,0.2,0.2))

tree<-read.nexus("Test_tree.tre")
crown.age<- max(node.age(tree)$ages)

###################
# a) Paleo-temprature
###################

InfTemp<-read.table("./PaleoEnv/PastTemperature.txt", header=T) 

if (crown.age>=100){max.y<-60}
if (crown.age<=100){max.y<-50}
if (crown.age<=66){max.y=22}
if (crown.age<=35){max.y=15}
if (crown.age<=20){max.y=15}

plot(-InfTemp[,"Age"],InfTemp[,"Temperature"],xlab="", ylab="Temperature (°C)", type="p", lwd=0.1, xlim=c(-crown.age,0), ylim=c(0, max.y), col="grey", las=1, cex=0.8, cex.axis=0.8, bty = "n")
temp.spl<-sm.spline(-InfTemp[,"Age"], InfTemp[,"Temperature"],df=100)
lines(temp.spl, col="firebrick1", lwd=2)
legend("topleft", bty="n", c("a) Paleoclimate"), cex=0.7)

abline(v=c(-298.9,-252.2,-201.3,-145,-100.5,-66,-56,-47.8,-33.9,-28.1,-23.03,-15.97,-11.62,-5.33,-2.58),col="grey",lty="dotted",lwd="1") # add vertical lines to delineate geological periods/epochs/stages


###################
# b) Paleo-altitude
###################

InfAlti<-read.table("./PaleoEnv/PastAndeanAltitude.txt", header=T)

plot(-InfAlti[,"Age"],InfAlti[,"Altitude"],xlab="", ylab="Altitude (m)", type="p", lwd=0.5, xlim=c(-crown.age,0), col="grey", las=1, cex=0.8, cex.axis=0.8, bty = "n")
alti.spl<-sm.spline(-InfAlti[,"Age"], InfAlti[,"Altitude"],df=20)
lines(alti.spl, col="burlywood4", lwd=2)
legend("topleft", bty="n", c("b) Andean orogeny"), cex=0.7)

abline(v=c(-66,-56,-47.8,-33.9,-28.1,-23.03,-15.97,-11.62,-5.33,-2.58),col="grey",lty="dotted",lwd="1") # add vertical lines to delineate geological periods/epochs/stages


###################
# c) Paleo-sea level
###################

InfSea<-read.table("./PaleoEnv/PastSeaLevel.txt", header=T) 

if (crown.age>=35){max.y=150}
if (crown.age>=35){min.y=-100}
if (crown.age<=35){max.y=50}
if (crown.age<=35){min.y=-100}

plot(-InfSea[,"Age"], InfSea[,"Level"],xlab="", ylab="Sea level (m)", type="p", lwd=0.1, xlim=c(-crown.age,0), ylim=c(min.y, max.y), col="grey", las=1, cex=0.8, cex.axis=0.8, bty = "n")
sea.spl<-sm.spline(-InfSea[,"Age"], InfSea[,"Level"],df=100)
lines(sea.spl, col="dodgerblue", lwd=2)
legend("topleft", bty="n", c("c) Past sea-level fluctuations"), cex=0.7)

abline(v=c(-298.9,-252.2,-201.3,-145,-100.5,-66,-56,-47.8,-33.9,-28.1,-23.03,-15.97,-11.62,-5.33,-2.58),col="grey",lty="dotted",lwd="1") # add vertical lines to delineate geological periods/epochs/stages


######################################
# d) Diversification as a function of temperature
######################################

attach("Test_treeM_complete_results_PastTemperature.Rdata")
mean.results<-final_table_tree_file[[1]][order(as.numeric(final_table_tree_file[[1]][,"AICc"])),]
std.err.results<-final_table_tree_file[[2]][order(as.numeric(final_table_tree_file[[2]][,"AICc"])),]

res<-sm.spline(InfTemp[,"Age"],InfTemp[,"Temperature"],df=100) # You may adjust the df value (the higher you put df the more accurate the curve will be)
Temp_fun<-function(x){predict(res,x)}

# Speciation rate
# Please replace the as.numeric() in the objects TempDep_lamb_par1 / 2 below with the parameter estimates of your paleoenvironment-dependent analysis. 
# TempDep_lamb_par1 is lambda, TempDep_lamb_par2 is alpha.
# TempDep_lamb_par1_sd is the standard error of lambda, and TempDep_lamb_par2_sd is the standard error of alpha. 

TempDep_lamb_par1 <- as.numeric(mean.results[1,5])
TempDep_lamb_par2 <- as.numeric(mean.results[1,6]); if (TempDep_lamb_par2 =="NaN"){TempDep_lamb_par2=0}
TempDep_lamb_par1_sd <- 0 # can be changed if you ran the analyses on a set of trees (then you have the standard errors table
TempDep_lamb_par2_sd <- 0 # can be changed if you ran the analyses on a set of trees (then you have the standard errors table)

f.lamb.mean<-function(x){TempDep_lamb_par1*exp(TempDep_lamb_par2*Temp_fun(x))}
f.lamb.low<-function(x){(TempDep_lamb_par1-TempDep_lamb_par1_sd)*exp((TempDep_lamb_par2-TempDep_lamb_par2_sd)*Temp_fun(x))}
f.lamb.high<-function(x){(TempDep_lamb_par1+TempDep_lamb_par1_sd)*exp((TempDep_lamb_par2+TempDep_lamb_par2_sd)*Temp_fun(x))}

max.y<-round(max(f.lamb.high(InfTemp[,"Age"])),1)

plot(-InfTemp[,"Age"], f.lamb.mean(InfTemp[,"Age"]), ty="l",col="chartreuse3",xlim=c(-crown.age,0), ylim=c(0,max.y), lwd=2, yaxt="n", xlab="Time (Myrs ago)",ylab="Speciation (green) and extinction (red) rates", cex.axis=0.8, bty = "n")
axis(2, at = seq(0, max.y, by = 0.05), las=1, cex.axis=0.8)
legend("topleft", bty="n", c("d) Diversification according to temperature"), cex=0.7)
lines(-InfTemp[,"Age"], f.lamb.low(InfTemp[,"Age"]),ty="l",col="chartreuse3",xlim=c(-crown.age,0),lwd=1,lty="dotted",yaxt="n")
lines(-InfTemp[,"Age"], f.lamb.high(InfTemp[,"Age"]),ty="l",col="chartreuse3",xlim=c(-crown.age,0),lwd=1,lty="dotted",yaxt="n")

# Extinction rate
# Please replace the as.numeric() in the objects TempDep_mu_par1 / 2 below with the parameter estimates of your paleoenvironment-dependent analysis. 
# TempDep_mu_par1 is mu, and TempDep_mu_par2 is beta.
# TempDep_mu_par1_sd is the standard error of mu, and TempDep_mu_par2_sd is the standard error of beta. 

TempDep_mu_par1 <- as.numeric(mean.results[1,7]); if (TempDep_mu_par1=="NaN"){TempDep_mu_par1=0}
TempDep_mu_par2 <- as.numeric(mean.results[1,8]); if (TempDep_mu_par2=="NaN"){TempDep_mu_par2=0}
TempDep_mu_par1_sd <- 0 # can be changed if you ran the analyses on a set of trees (then you have the standard errors table
TempDep_mu_par2_sd <- 0 # can be changed if you ran the analyses on a set of trees (then you have the standard errors table

f.mu.mean<-function(x){TempDep_mu_par1*exp(TempDep_mu_par2*Temp_fun(x))}
f.mu.low<-function(x){(TempDep_mu_par1-TempDep_mu_par1_sd)*exp((TempDep_mu_par2-TempDep_mu_par2_sd)*Temp_fun(x))}
f.mu.high<-function(x){(TempDep_mu_par1 + TempDep_mu_par1_sd)*exp((TempDep_mu_par2 + TempDep_mu_par2_sd)*Temp_fun(x))}

lines(-InfTemp[,"Age"], f.mu.mean(InfTemp[,"Age"]), ty="l",col="red",lwd=1, yaxt="n", xlab="Time (Myrs ago)",ylab="Speciation (blue) and extinction (red) rates", cex.axis=0.8, bty = "n")
lines(-InfTemp[,"Age"], f.mu.low(InfTemp[,"Age"]),ty="l",col="red",lwd=1,lty="dotted",yaxt="n")
lines(-InfTemp[,"Age"], f.mu.high(InfTemp[,"Age"]),ty="l",col="red",lwd=1,lty="dotted",yaxt="n")

abline(v=c(-298.9,-252.2,-201.3,-145,-100.5,-66,-56,-47.8,-33.9,-28.1,-23.03,-15.97,-11.62,-5.33,-2.58),col="grey",lty="dotted",lwd="1") # add vertical lines to delineate geological periods/epochs/stages


######################################
# e) Diversification as a function of Andean altitude
######################################

attach("Test_treeM_complete_results_PastAndeanAltitude.Rdata")
mean.results<-final_table_tree_file[[1]][order(as.numeric(final_table_tree_file[[1]][,"AICc"])),]
std.err.results<-final_table_tree_file[[2]][order(as.numeric(final_table_tree_file[[2]][,"AICc"])),]

res<-sm.spline(InfAlti[,"Age"], InfAlti[,"Altitude"],df=20)
Alti_fun <-function(x){predict(res,x)}

# Please replace the as.numeric() in the objects AltiDep_lamb_par1 / 2 below with the parameter estimates of your paleoenvironment-dependent analysis. 
# TempDep_lamb_par1 is lambda, TempDep_lamb_par2 is alpha.
# TempDep_lamb_par1_sd is the standard error of lambda, and TempDep_lamb_par2_sd is the standard error of alpha. 

# Speciation is varying:
AltiDep_lamb_par1 <- as.numeric(mean.results[1,5])
AltiDep_lamb_par2 <- as.numeric(mean.results[1,6]); if (TempDep_lamb_par2 =="NaN"){TempDep_lamb_par2=0}
AltiDep_lamb_par1_sd <- 0
AltiDep_lamb_par2_sd <- 0

f.lamb.mean<-function(x){AltiDep_lamb_par1*exp(AltiDep_lamb_par2*Alti_fun(x))}
f.lamb.low<-function(x){(AltiDep_lamb_par1-AltiDep_lamb_par1_sd)*exp((AltiDep_lamb_par2-AltiDep_lamb_par2_sd)*Alti_fun(x))}
f.lamb.high<-function(x){(AltiDep_lamb_par1+AltiDep_lamb_par1_sd)*exp((AltiDep_lamb_par2+AltiDep_lamb_par2_sd)*Alti_fun(x))}

max.y<-round(max(f.lamb.high(InfAlti[,"Age"])),1)

plot(-InfAlti[,"Age"], f.lamb.mean(InfAlti[,"Age"]), ty="l",col="chartreuse3",xlim=c(-crown.age,0), ylim=c(0,max.y), lwd=2, yaxt="n", xlab="Time (Myrs ago)",ylab="Speciation (green) and extinction (red) rates", cex.axis=0.8, bty = "n")
axis(2, at = seq(0, max.y, by = 0.05), las=1, cex.axis=0.8)
legend("topleft", bty="n", c("e) Diversification according to Andean altitude"), cex=0.7)
lines(-InfAlti[,"Age"], f.lamb.low(InfAlti[,"Age"]),ty="l",col="chartreuse3",lwd=1,lty="dotted",yaxt="n")
lines(-InfAlti[,"Age"], f.lamb.high(InfAlti[,"Age"]),ty="l",col="chartreuse3",lwd=1,lty="dotted",yaxt="n")


# Extinction rate
# Please replace the as.numeric() in the objects AltiDep_mu_par1 / 2 below with the parameter estimates of your paleoenvironment-dependent analysis. 
# TempDep_mu_par1 is mu, and TempDep_mu_par2 is beta.
# TempDep_mu_par1_sd is the standard error of mu, and TempDep_mu_par2_sd is the standard error of beta. 

AltiDep_mu_par1 <- as.numeric(mean.results[1,7]); if (TempDep_mu_par1=="NaN"){TempDep_mu_par1=0}
AltiDep_mu_par2 <- as.numeric(mean.results[1,7]); if (TempDep_mu_par2=="NaN"){TempDep_mu_par2=0}
AltiDep_mu_par1_sd <- 0
AltiDep_mu_par2_sd <- 0

f.mu.mean<-function(x){AltiDep_mu_par1*exp(AltiDep_mu_par2*Alti_fun(x))}
f.mu.low<-function(x){(AltiDep_mu_par1-AltiDep_mu_par1_sd)*exp((AltiDep_mu_par2-AltiDep_mu_par2_sd)*Alti_fun(x))}
f.mu.high<-function(x){(AltiDep_mu_par1+AltiDep_mu_par1_sd)*exp((AltiDep_mu_par2+AltiDep_mu_par2_sd)*Alti_fun(x))}

lines(-InfAlti[,"Age"], f.mu.mean(InfAlti[,"Age"]), ty="l",col="red",xlim=c(-crown.age,0), ylim=c(0,0.8), lwd=2, yaxt="n", xlab="Time (Myrs ago)",ylab="", cex.axis=0.8, bty = "n")
lines(-InfAlti[,"Age"], f.mu.low(InfAlti[,"Age"]),ty="l",col="red",lwd=1,lty="dotted",yaxt="n")
lines(-InfAlti[,"Age"], f.mu.high(InfAlti[,"Age"]),ty="l",col="red",lwd=1,lty="dotted",yaxt="n")

abline(v=c(-66,-56,-47.8,-33.9,-28.1,-23.03,-15.97,-11.62,-5.33,-2.58),col="grey",lty="dotted",lwd="1") # add vertical lines to delineate geological periods/epochs/stages


######################################
# f) Diversification as a function of sea level
######################################

attach("Test_treeM_complete_results_PastSeaLevel.Rdata")
mean.results<-final_table_tree_file[[1]][order(as.numeric(final_table_tree_file[[1]][,"AICc"])),]
std.err.results<-final_table_tree_file[[2]][order(as.numeric(final_table_tree_file[[2]][,"AICc"])),]

res<-sm.spline(InfSea[,"Age"], InfSea[,"Level"],df=50)
Sea_fun <-function(x){predict(res,x)}

# Please replace the as.numeric() in the objects SeaDep_lamb_par1 / 2 below with the parameter estimates of your paleoenvironment-dependent analysis. 
# SeaDep_lamb_par1 is lambda, SeaDep_lamb_par2 is alpha.
# SeaDep_lamb_par1_sd is the standard error of lambda, and SeaDep_lamb_par2_sd is the standard error of alpha. 

# Speciation is varying:
SeaDep_lamb_par1 <- as.numeric(mean.results[1,5])
SeaDep_lamb_par2 <- as.numeric(mean.results[1,6]); if (TempDep_lamb_par2 =="NaN"){TempDep_lamb_par2=0}
SeaDep_lamb_par1_sd <- 0
SeaDep_lamb_par2_sd <- 0

f.lamb.mean<-function(x){SeaDep_lamb_par1*exp(SeaDep_lamb_par2*Sea_fun(x))}
f.lamb.low<-function(x){(SeaDep_lamb_par1-SeaDep_lamb_par1_sd)*exp((SeaDep_lamb_par2-SeaDep_lamb_par2_sd)*Sea_fun(x))}
f.lamb.high<-function(x){(SeaDep_lamb_par1+SeaDep_lamb_par1_sd)*exp((SeaDep_lamb_par2+SeaDep_lamb_par2_sd)*Sea_fun(x))}

max.y<-round(max(f.lamb.high(InfSea[,"Age"])),1)

plot(-InfSea[,"Age"], f.lamb.mean(InfSea[,"Age"]), ty="l",col="chartreuse3",xlim=c(-crown.age,0), ylim=c(0,max.y), lwd=2, yaxt="n", xlab="Time (Myrs ago)",ylab="Speciation (green) and extinction (red) rates", cex.axis=0.8, bty = "n")
axis(2, at = seq(0, max.y, by = 0.05), las=1, cex.axis=0.8)
legend("topleft", bty="n", c("f) Diversification according to sea level"), cex=0.7)
lines(-InfSea[,"Age"], f.lamb.low(InfSea[,"Age"]),ty="l",col="chartreuse3",lwd=1,lty="dotted",yaxt="n")
lines(-InfSea[,"Age"], f.lamb.high(InfSea[,"Age"]),ty="l",col="chartreuse3",lwd=1,lty="dotted",yaxt="n")

# Extinction rate
# Please replace the as.numeric() in the objects SeaDep_mu_par1 / 2 below with the parameter estimates of your paleoenvironment-dependent analysis. 
# SeaDep_mu_par1 is mu, and SeaDep_mu_par2 is beta.
# SeaDep_mu_par1_sd is the standard error of mu, and SeaDep_mu_par2_sd is the standard error of beta. 

SeaDep_mu_par1 <- as.numeric(mean.results[1,7]); if (TempDep_mu_par1=="NaN"){TempDep_mu_par1=0}
SeaDep_mu_par2 <- as.numeric(mean.results[1,7]); if (TempDep_mu_par2=="NaN"){TempDep_mu_par2=0}
SeaDep_mu_par1_sd <- 0
SeaDep_mu_par2_sd <- 0

f.mu.mean<-function(x){SeaDep_mu_par1*exp(SeaDep_mu_par2*Alti_fun(x))}
f.mu.low<-function(x){(SeaDep_mu_par1-SeaDep_mu_par1_sd)*exp((SeaDep_mu_par2-SeaDep_mu_par2_sd)*Alti_fun(x))}
f.mu.high<-function(x){(SeaDep_mu_par1+SeaDep_mu_par1_sd)*exp((SeaDep_mu_par2+SeaDep_mu_par2_sd)*Alti_fun(x))}

lines(-InfSea[,"Age"], f.mu.mean(InfSea[,"Age"]), ty="l",col="red",xlim=c(-crown.age,0), lwd=2, yaxt="n", xlab="Time (Myrs ago)",ylab="", cex.axis=0.8, bty = "n")
lines(-InfSea[,"Age"], f.mu.low(InfSea[,"Age"]),ty="l",col="red",lwd=1,lty="dotted",yaxt="n")
lines(-InfSea[,"Age"], f.mu.high(InfSea[,"Age"]),ty="l",col="red",lwd=1,lty="dotted",yaxt="n")

abline(v=c(-66,-56,-47.8,-33.9,-28.1,-23.03,-15.97,-11.62,-5.33,-2.58),col="grey",lty="dotted",lwd="1") # add vertical lines to delineate geological periods/epochs/stages

dev.off()
