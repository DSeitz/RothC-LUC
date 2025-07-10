### RothC-LUC function with additional AGG pool (here: aggr) 
# to model SOC dynamics after land-use change between cropland and grassland as described in Seitz et al. 2025 (EJSS)
# The AGG pool has a maximum capacity (= aggr.sat) calculated via the difference between site-specific equilibrium SOC stocks: (SOC_grass - SOC_crop)*0.8
# The decomposition rate of the AGG pool (= aggr.fast) was set to 0.05
# The pool transfer of the AGG pool (= aggr.slow) was set to 0.3 in the case of grassland, and set to 0 in the case of cropland (line 135)
# default land use (= landuse) is cropland, making the model work as the default RothC
# packages needed: 
# library(SoilR)
# library(Rcpp)

sourceCpp(paste0(root.path,"/RothC_LUC.cpp")) # enter your path (root.path) here to use the cpp function

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
RothCaggrSoilRModel=function (t, res= 12,
                              ks = c(k.DPM = 10, k.RPM = 0.3,k.aggr=0.3, k.BIO = 0.66, k.HUM = 0.02, k.IOM = 0),
                              C0 = c(1, 2, 1,5, 111, 2.7),
                              In = c(1,1,0,0,0,0),
                              clay = 23.4, xi = 1,
                              aggr.sat, # maximum capacity of the aggr pool: (SOC_grass - SOC_crop)*0.8
                              aggr.fast, # decomposition rate: 0.05
                              aggr.slow, # pool transfer rate: 0.3
                              landuse="cropland",
                              solver = deSolve.lsoda.wrapper, pass = FALSE) 
{
  t.old=t
  t=seq.default(from = t[1],to = t[length(t)],by = 1/res)
  
  t_start = min(t)
  t_end = max(t)
  if (length(ks) != 6) 
    stop("ks must be of length = 6")
  if (length(C0) != 6) 
    stop("the vector with initial conditions must be of length = 5")
  ####################################################################
  # conditions for In
  ####################################################################
  if (class(In) == "data.frame") 
  {
    if(length(In)!=6)
    {
      stop("In RothC_LUCModel input must be exactly a dataframe with 6 cols")
    }
    if(length(In)!=6)
    {
      stop("In RothC_LUCModel input must be exactly a dataframe with 6 cols and the length of t in monthly resolution")
    }
    if(length(In[,1])!=length(t))
    {
      stop("In RothC_LUCModel input must be exactly a dataframe with 6 cols and the length of t in monthly resolution")
    }
  }
  if (class(In) == "numeric") 
  {
    if(length(In)!=6)
    {
      stop("In RothC_LUCModel input must be exactly a vector with 6 annual inputs (or a data.frame with 6 cols)")
    }
    dpm.in=rep(In[1]/res,length(t))
    rpm.in=rep(In[2]/res,length(t))
    aggr.in=rep(In[3]/res,length(t))
    bio.in=rep(In[4]/res,length(t))
    hum.in=rep(In[5]/res,length(t))
    iom.in=rep(In[6]/res,length(t))
    In=data.frame(dpm.in, rpm.in,aggr.in, bio.in,hum.in,iom.in)
  }
  
  if (class(xi)=="numeric")
  {
    if(length(xi)!=1)
    {
      stop("RothC_LUC mean responses must be exactely one value")
    }
    xi=rep(xi,length(t))
  }
  if (class(xi)=="data.frame")
  {
    if(length(xi)!=2)
    {
      stop("RothC_LUC: if xi is a data.frame it has exactly two cols first time second values")
    }
    if (sum(xi[,1]!=t)>0)
    {
      stop("RothC_LUC first column of xi -the time- is strange.")
    }
    xi=xi[,2]
  }
  #######################################################################################################################
  # here the transfer between aggregated and not aggregated stuff
  #######################################################################################################################
  c.pools=dataframe(ncol=6,nrow=length(t),fill = 0)
  names(c.pools)=c("c_dpm","c_rpm","c_aggr","c_bio","c_hum","c_iom")
  c.pools[1,]=C0
  
  
  x = 1.67 * (1.85 + 1.6 * exp(-0.0786 * clay))
  noco2 = 1-x/(x+1)
  ks.new = unlist(unname(ks))
  
  out = as.data.frame(RothCagg2_sub_cpp(t = t,ks = ks.new,xi = xi,In = as.matrix(In),c_pools = as.matrix(c.pools),
                          noco2 = noco2,res = res,aggr_sat = aggr.sat,aggr_fast = aggr.fast,
                          aggr_slow = aggr.slow,landuse = landuse))
 
  out1=out[t %in% t.old,]
  out1
}


