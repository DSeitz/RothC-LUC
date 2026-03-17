# RothC-LUC to model land-use change
# Copyright (C) 2026
# by
# René Dechow, Thünen Institute of Climate-Smart Agriculture, Braunschweig, Germany
# Daria Seitz, Thünen Institute of Climate-Smart Agriculture, Braunschweig, Germany

### Example script to illustrate differences in the pool dynamics between RothC-default and RothC-LUC
# RothC-LUC to model land-use change effects on soil organic carbon

rm(list = ls())
library(Rcpp)
library(SoilR)
root.path = "/path-to-folder/RothC-LUC-main/" # enter your path (root.path) here to use the cpp function
source(paste0(root.path,"RothC_LUC_Model.R"))
sourceCpp(paste0(root.path,"/RothC_LUC.cpp")) 


func.plot.timeseries=function(t,out.default,out.aggr)
{
  par(mfrow = c(1,1))
  colors = c("red","orange","blue","green","brown","grey","black")
  matplot(t,out.default,col = colors[c(1,2,4,5,6,7)], lwd =2, type ="l",lty=1,ylim= c(-20,50),
          ylab =expression("Corg [Mg ha"^"-1"*"]"))
  matplot(t,out.aggr,type = "l",col= colors,lty=2,add = TRUE)
  legend(title = "RothC-default",x = "bottomleft",legend = c("dpm","rpm","bio","hum","iom","Ctot"),col =colors[c(1,2,4,5,6,7)],bty ="n",ncol=2,lty=1)
  legend(title = "RothC-LUC",x = "bottomright",legend = c("dpm","rpm","aggr","bio","hum","iom","Ctot"),col =colors,bty ="n",ncol=3,lty=2)
  
}


################################################################################
# define model input and parameters
################################################################################
aggr.slow = 0.3 # fraction of decomposition products entering the AGGR pool, if landuse == grassland
aggr.fast = 0.05 # decomposition rate of the AGGR pool
aggr.sat = 11 # capacity of AGGR pool [t C/ha]


################################################################################
# general settings
################################################################################
t = seq.default(from = 0,to = 100,by = 1/12)# hundred years
c_input = 3 # t C / ha
DR = 3
xi = 0.6 # rate modifying factor
clay = 10 #clay content
################################################################################
# settings for SoilR RothC
################################################################################
C0.ori =c(0,3,0.2,32,3) #initial conditions
################################################################################
# Settings RothC aggr
################################################################################
In = c(c_input*(DR/(DR + 1)),c_input * (1/(DR + 1)),0,0,0,0) # annual inputs into pool
C0 =c(0,3,3,0.2,32,3) #initial conditions

################################################################################
# run RothC_SoilR 
################################################################################
rothc.default = RothCModel(t = t,C0 = C0.ori,xi = xi,In = c_input,clay = clay,DR = DR)
out.default = getC(rothc.default)
out.default = data.frame(out.default)
out.default = data.frame(out.default, rowSums(out.default))


################################################################################
# 1. run the example assuming cropland
################################################################################
landuse = "cropland"

out.aggr =RothCAggrRModel(t=t,In = In,C0 = C0,xi = xi,aggr.fast = aggr.fast,aggr.slow = aggr.slow, aggr.sat = aggr.sat,landuse = landuse,clay = clay )
out.aggr = data.frame(out.aggr)
out.aggr = data.frame(out.aggr, rowSums(out.aggr))

func.plot.timeseries(t = t,out.default = out.default,out.aggr = out.aggr)
################################################################################
# 2. run the example assuming grassland
################################################################################
landuse = "grassland"

out.aggr =RothCAggrRModel(t=t,In = In,C0 = C0,xi = xi,aggr.fast = aggr.fast,aggr.slow = aggr.slow, aggr.sat = aggr.sat,landuse = landuse,clay = clay )
out.aggr = data.frame(out.aggr)
out.aggr = data.frame(out.aggr, rowSums(out.aggr))

func.plot.timeseries(t = t,out.default = out.default,out.aggr = out.aggr)
################################################################################
# 3. run the example assuming landuse change cropland - grassland after 20 years
################################################################################
landuse = rep("grassland",length(t))
landuse[1:240] = "cropland"

out.aggr =RothCAggrRModel(t=t,In = In,C0 = C0,xi = xi,aggr.fast = aggr.fast,aggr.slow = aggr.slow, aggr.sat = aggr.sat,landuse = landuse,clay = clay )
out.aggr = data.frame(out.aggr)
out.aggr = data.frame(out.aggr, rowSums(out.aggr))

func.plot.timeseries(t = t,out.default = out.default,out.aggr = out.aggr)
################################################################################
# 4. run the example assuming landuse change grassland - cropland after 20 years
################################################################################
landuse = rep("cropland",length(t))
landuse[1:240] = "grassland"

out.aggr =RothCAggrRModel(t=t,In = In,C0 = C0,xi = xi,aggr.fast = aggr.fast,aggr.slow = aggr.slow, aggr.sat = aggr.sat,landuse = landuse,clay = clay )
out.aggr = data.frame(out.aggr)
out.aggr = data.frame(out.aggr, rowSums(out.aggr))

func.plot.timeseries(t = t,out.default = out.default,out.aggr = out.aggr)

