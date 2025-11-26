###############################################################################
# Example of model initialisation in R : equilibrium run
# coded and provided by René Dechow (Thünen Institute)
# There  are two model runs, the first model is initialized by an equilibrium run
# The second model run just by fixed values


rm(list=ls())
library(SoilR)

root.base = "D:/folder_with_script/" ## enter your path here
setwd(root.base)
root=paste0(root.base,"subfolder_with_data/") ## enter your path here
setwd(root)
cols=c("red","orange","green","brown","grey","black")
font.scale=1.5
################################################################################
# set parameters 
################################################################################
#  Define Tillage options
options = c("conventional","reduced","no")
rmf.Till =c( 1,             0.93, 0.95)
till.cats= data.frame (options,rmf.Till)

################################################################################
# Definition of the simulation period and depth
################################################################################
years=seq(1/12,500,by=1/12) # simulation period in monthly time steps
soil.thick=25  #Soil thickness (topsoil), in cm
################################################################################
# site properties
################################################################################
SOC=69.7       #Soil organic carbon in Mg/ha 
clay=48        #Percent clay
################################################################################
# Definition of soil management
################################################################################
Cinputs1 = 2.7   #Annual C inputs to soil in Mg/ha/yr
Cinputs2 = 3   #Annual C inputs to soil in Mg/ha/yr

tillage1 = "conventional"
tillage2 = "reduced"
 
bare1 = TRUE
bare2 = TRUE
################################################################################
# get the weather data
################################################################################
file.name = "weather-data-example.txt"
clim = read.table(file = file.name,header = TRUE,sep = "\t")


################################################################################
# Quantification of "Rate Modulating Factors"
################################################################################
rmf.T = fT.RothC(clim$Temp) #Temperature effects per month
rmf.W1 = fW.RothC(P = clim$Precip, 
            E= clim$Evp, 
            S.Thick = soil.thick, 
            pClay = clay, 
            pE = 1.0, bare = bare1)$b #Moisture effects per month
rmf.Cover = rep(ifelse (bare1,1,0.6),length(rmf.T))
rmf.Till1  = rep(till.cats$rmf.Till[till.cats$options==tillage1],length(rmf.T))


rmf.all1 = rmf.T*rmf.W1*rmf.Cover*rmf.Till
xi.frame1 =data.frame(years,rep(rmf.all1,length.out=length(years)))

rmf.T = fT.RothC(clim$Temp) #Temperature effects per month
rmf.W2 = fW.RothC(P = clim$Precip, 
                 E= clim$Evp, 
                 S.Thick = soil.thick, 
                 pClay = clay, 
                 pE = 1.0, bare = bare1)$b #Moisture effects per month
rmf.Cover = rep(ifelse (bare1,1,0.6),length(rmf.T))
rmf.Till2  = rep(till.cats$rmf.Till[till.cats$options==tillage2],length(rmf.T))


rmf.all2 = rmf.T*rmf.W2*rmf.Cover*rmf.Till
xi.frame2 =data.frame(years,rep(rmf.all2,length.out=length(years)))

################################################################################
# define assumption for the initial pool distribution
FallIOM =0.049*50^(1.139) #IOM using Falloon method
C0_1=c(DPM=0, RPM= 50, BIO=0, HUM=0, IOM=FallIOM)
C0_2=c(DPM= 0, RPM= 50, BIO=0, HUM=0, IOM=FallIOM) 
################################################################################
# Model initialisation function
################################################################################
fit_the_inital_C = function (In, C0 , clay, xi , C.obs, calli.mod)
{
  Model_init=RothCModel(t= 0:3000,C0 = C0_1,
                    In=In, clay=clay, xi=xi) #Loads the model
  Ct.in = getC(Model_init) #Calculates stocks for each pool per month
  C.mod = sum(Ct.in[nrow(Ct.in),])
  
  output.list =list()
  
  fit = (C.obs-C.mod)^2
  init.C = Ct.in[nrow(Ct.in),]
  
  #output.list
  if (calli.mod == 1)
  {
    retour = fit
  }else{
    retour = init.C
  }
  retour
}

#################################
# Parameter for initialisation
#################################
cin = 2 # start value
cin.min = 0.1 # lowest value possible (assumption)
cin.max = 20  # highest value possible (assumption)
xi.cin = mean (xi.frame1[,2]) # assumed boundary conditions
C.obs = 90
calli.mod = 1

################################################################################
# Optimizing the initial pool distribution
################################################################################
fit.par = optim(par = cin,fn = fit_the_inital_C, cin, C0_1, clay, xi.cin, C.obs,calli.mod,
                method = "Brent",lower = cin.min,upper = cin.max)
C0_equi = C0_1
C0_equi = fit_the_inital_C (In = fit.par$par, C0 = C0_1, clay = clay, xi = xi.cin , C.obs = C.obs,calli.mod = 0)


################################################################################
# Model application
################################################################################
Model1=RothCModel(t=years,C0 = C0_equi,
                  In=Cinputs1, clay=clay, xi=xi.frame1) #Loads the model
Ct1=getC(Model1) #Calculates stocks for each pool per month
Ct1 = data.frame(Ct1)
Ctot = rowSums(Ct1)
Ct1 = data.frame(Ct1,Ctot)

Model2=RothCModel(t=years,C0 = C0_2,
                  In=Cinputs1, clay=clay, xi=xi.frame1) #Loads the model
Ct2=getC(Model2) #Calculates stocks for each pool per month
Ct2 = data.frame(Ct2)
Ctot = rowSums(Ct2)
Ct2 = data.frame(Ct2,Ctot)
################################################################################
# plot results
################################################################################

par(mfrow=c(1,1))
matplot(years, Ct1, type="l", lty=1, col= cols,ylim=c(0,190),lwd=3,
        xlab="Time (years)", ylab="C stocks (Mg/ha)",cex.axis=font.scale,cex.lab=font.scale)
matplot(years, Ct2, type="l", lty=2, col= cols,add=TRUE,lwd=0.5,
        xlab="Time (years)", ylab="C stocks (Mg/ha)",)

legend("topleft", c("DPM", "RPM", "BIO", "HUM", "IOM","Corg"),
       lty=1, col= cols, bty="n",cex=font.scale,ncol = 2)
# What you can see here:
# Ct1 is initialized by an equilibrium run (fat line) (most C was allocated to the HUM pool during initialization)
# Ct2 is run just by fixed values (with most C in the RPM pool)

