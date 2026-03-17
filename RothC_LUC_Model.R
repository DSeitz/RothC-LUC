# RothC-LUC to model land-use change
# Copyright (C) 2026
# by
# René Dechow, Thünen Institute of Climate-Smart Agriculture, Braunschweig, Germany
# Daria Seitz, Thünen Institute of Climate-Smart Agriculture, Braunschweig, Germany


### RothC-LUC function with additional AGG pool (here: aggr) 
# to model SOC dynamics after land-use change between cropland and grassland as described in Seitz et al. 2025 (EJSS)
# The AGG pool has a maximum capacity (= aggr.sat) calculated via the difference between site-specific equilibrium SOC stocks: (SOC_grass - SOC_crop)*0.8
# The decomposition rate of the AGG pool (= aggr.fast) was set to 0.05
# The pool transfer of the AGG pool (= aggr.slow) was set to 0.3 in the case of grassland, and set to 0 in the case of cropland (line 135)
# default land use (= landuse) is cropland, making the model work as the default RothC
# packages needed: 
# library(SoilR)
# library(Rcpp)

# function: create a dataframe
dataframe<-function(ncol=10,nrow=10,basename="basename",fill=-9999)
{
  a=rep(fill,nrow)
  b=data.frame(a)
  if(ncol>1)
  {
    for(i in 2:ncol)
    {
      b=cbind(b,a)
    }
  }
  nam=rep(basename,ncol)
  nums=1:ncol
  nam=paste(nam,nums,sep="")
  names(b)=nam
  b
}

# function: RothC-LUC
RothCAggrRModel=function (t, ks= c(k.DPM = 10, k.RPM = 0.3,k.aggr = aggr.fast, k.BIO = 0.66, k.HUM = 0.02,k.IOM = 0), C0 = c(0, 0, 0,0, 0, 2.7), In=1.7,xi=1,
                          clay = 23.4,aggr.fast=aggr.fast,aggr.slow= aggr.slow,aggr.sat = aggr.sat,landuse ="cropland") 
{
  n.pool = length(ks)
  ks[3] = aggr.fast
  if (length(t)<3)
  {
    tseq =seq.default(from = t[1],to = t[2],by = 1/12)
  }else{
    tseq = t
  }
  if (class(In) ==  "numeric")
  {
    In.frame = dataframe(ncol = n.pool, nrow = length(tseq))
    if (length(In) != n.pool) stop("Length of In must correspond to N of pools.\n")
    for (i in 1:length(In))
    {
      In.frame[,i] = In[i]/12
    }
    In = In.frame
  } else if (class(In) == "data.frame"){
    if (ncol(In) != n.pool) stop ("Col.N of In do not correspond  to number of pools\n ")
    if (nroe(In) != n.pool) stop ("Row.N of In do not correspond time steps\n ")
  }
  
  if ((class(xi) == "numeric")& (length(xi) == 1))
  {
    xi.frame = dataframe(ncol = n.pool,nrow = length(tseq),fill = xi)
    xi.frame[,3] = 1
    xi = xi.frame
  }else if((class(xi) == "numeric")& (length(xi) == n.pool)){
    xi.frame = dataframe(ncol = n.pool,nrow = length(tseq))
    for (i in 1:length(xi.frame))
    {
      xi.frame[,i] = xi[i]
    }
    xi = xi.frame
    xi[,3] = 1
  }else if ((class(xi) == "data.frame") & (ncol(xi) == n.pool) & nrow(xi) == length(tseq)){
    xi[,3] = 1
  }else{
    stop("Error in format of xi\n")
  }
  #############################
  # landuse vector
  ############################
  if (length(landuse) == 1)
  {
    is.crop = rep(landuse =="cropland",length(tseq))
  }else{
    if (length(landuse) != length(tseq)) stop("Landuse vector length does not match N of time steps\n")
    is.crop = (landuse == "cropland")
  }
  
  
  ###############################################################
  # build the matrix grassland
  ###############################################################
  aggr.rest = 1-aggr.slow
  x = 1.67 * (1.85 + 1.6 * exp(-0.0786 * clay))
  B = aggr.rest * 0.46/(x + 1)
  H = aggr.rest * 0.54/(x + 1)
  ai3 = aggr.slow * ks
  ai4 = B * ks
  ai5 = H * ks
  A.grass = diag(-ks)
  A.grass[4, ] = A.grass[4, ] + ai4
  A.grass[5, ] = A.grass[5, ] + ai5
  A.grass[3, ] = A.grass[3, ] + ai3
  ###############################################################
  # build the matrix cropland
  ###############################################################
  aggr.trans = 0
  aggr.rest = 1-0
  x = 1.67 * (1.85 + 1.6 * exp(-0.0786 * clay))
  B = aggr.rest * 0.46/(x + 1)
  H = aggr.rest *0.54/(x + 1)
  ai4 = B * ks
  ai5 = H * ks
  A.crop = diag(-ks)
  A.crop[4, ] = A.crop[4, ] + ai4
  A.crop[5, ] = A.crop[5, ] + ai5
  ###############################################################
  # set Cin matrix
  #################### ###########################################
  Cin=as.matrix(In[,1:length(In)])/1
  ###############################################################
  # build xii matrix
  ###############################################################
  xii=as.matrix(xi)
  ###############################################################
  # Run the model
  ###############################################################
  out= as.data.frame(rothc_aggr_cpp(tseq = tseq,C0 = C0,Cin = Cin,xi = xii,A_crop = A.crop ,A_grass=A.grass, aggr_sat =aggr.sat,is_crop = is.crop ))
  out1=out[tseq %in% t,]
  names(out1) = c("c_dpm","c_rpm","c_aggr","c_bio","c_hum","c_iom")
  out1
}



