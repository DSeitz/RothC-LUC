
// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix RothCagg2_sub_cpp(
  NumericVector t,
  NumericVector ks, 
  NumericVector xi,      
  NumericMatrix In,
  NumericMatrix c_pools,
  double noco2,
  double res,
  double aggr_sat,
  double aggr_fast,
  double aggr_slow,
  String landuse
) 
{

   NumericVector c_loss(6);
  NumericVector minops = {0,0};
  double transfer;
  double transfer_agg;
  double transfer_rate_agg;
  int nsteps=t.length();
  int i,no;
  // time loop
  for(no = 1; no<nsteps; no++) 
  {  
    // dpm
    c_pools(no,0)=In(no,0)+c_pools(no-1,0);
    //rpm
    c_pools(no,1)=In(no,1)+c_pools(no-1,1);
    //hum
    c_pools(no,3)=In(no,3)+c_pools(no-1,3);
    // all the rest
    c_pools(no,2)=In(no,2)+c_pools(no-1,2);
    
    c_pools(no,4)=In(no,4)+c_pools(no-1,4);
    c_pools(no,5)=In(no,5)+c_pools(no-1,5);
    //#####################################################
    //  # set the ks.aggr in dependence on landuse and aggr.sat
    //#####################################################
    if((landuse=="grassland") & (c_pools(no,3)<aggr_sat))
    {
      ks[2]=aggr_fast;
      transfer_rate_agg=aggr_slow;
    }else{
      ks[2]=aggr_fast;
      transfer_rate_agg=0;
    }
    //#####################################################
    //#bulk soil decomposition
    //#####################################################
    
    for (i=0;i<6;i++)
    {
      c_loss[i]=c_pools(no,i)*(1-exp(-ks[i]/res*xi[no]));
      c_pools(no,i)=c_pools(no,i)-c_loss[i];
    }
    minops[0] =  (c_loss[0]+c_loss[1]+c_loss[2]+c_loss[3]+c_loss[4])*transfer_rate_agg;
    minops[1] =  aggr_sat-c_pools(no,2);               
    transfer_agg = min(minops);
    transfer=((c_loss[0]+c_loss[1]+c_loss[2]+c_loss[3]+c_loss[4])-transfer_agg)*noco2;
    //################################################
    //# distribution of C  among pools BIO and HUM
    //################################################
    c_pools(no,2)=transfer_agg+c_pools(no,2);
    c_pools(no,3)=(0.46*transfer)+c_pools(no,3);
    c_pools(no,4)=(0.54*transfer)+c_pools(no,4);
  }
  return c_pools;
}
  

  