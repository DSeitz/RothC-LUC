  // RothC-LUC February 2026
  // by 
  // René Dechow, Thünen Institute of Climate-Smart Agriculture, Braunschweig, Germany
  // Daria Seitz, Thünen Institute of Climate-Smart Agriculture, Braunschweig, Germany
  //
  // modified from:
  // sorcering V.1.0 March 2021
  // soil organic carbon & nitrogen modelling groundwork 
  // by 
  // Dr. Marc Scherstjanoi, Thünen Institute of Forest Ecosystems, Eberswalde, Germany &
  // René Dechow, Thünen Institute of Climate-Smart Agriculture, Braunschweig, Germany



  
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
List
rothc_aggr_cpp(
  const Nullable<NumericVector&> tseq = R_NilValue,
  const bool tautomatic = false, 
  const Nullable<NumericVector&> C0 = R_NilValue, 
  const Nullable<LogicalVector&> is_crop = R_NilValue, 
  const Nullable<NumericMatrix&> Cin = R_NilValue,
  const Nullable<NumericMatrix&> xi = R_NilValue,
  const Nullable<NumericMatrix&> A_crop = R_NilValue, 
  const Nullable<NumericMatrix&> A_grass = R_NilValue, 
  double aggr_sat = 10,
  const bool meanCin = false 
) 
{    
  //  Environment base = Environment("package:base");
  //  Function readline = base["readline"];
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //input transformation
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  int timediv = 12;

  // transfer matrix 
  NumericMatrix Ac_(A_crop);
  NumericMatrix Ag_(A_grass); 
  
  arma::mat Ac_arma = as<arma::mat>(Ac_);
  arma::mat Ag_arma = as<arma::mat>(Ag_);
  arma::mat A_arma = Ac_arma;
  int nr_pools = Ac_arma.n_rows; // A defines nr_pools    
  
  // input, starting conditions, environment
  NumericVector C0_ (nr_pools); 
  if (C0.isNotNull()) C0_=C0;
  NumericVector tseq_  = {0,1};
  if (tseq.isNotNull()) tseq_=tseq;
  LogicalVector is_crop_ = {FALSE,FALSE};
  if (is_crop.isNotNull()) is_crop_ = is_crop;
  
  NumericMatrix Cin_(tseq_.length(),nr_pools); 
  NumericVector v (tseq_.length()*nr_pools,1);
  NumericMatrix xi_(tseq_.length(),nr_pools,v.begin());
  
  if (Cin.isNotNull()) 
  {
    NumericMatrix CinX(Cin);
    Cin_=CinX;
  }
  if (xi.isNotNull()) 
  {   
    NumericMatrix xiX(xi);
    xi_=xiX;
  }
  
  // Rccp to RccpArmadillo
  arma::vec C0_arma = as<arma::vec>(C0_);
  arma::mat Cin_arma = as<arma::mat>(Cin_);
  arma::mat xi_arma = as<arma::mat>(xi_);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //input verification
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  String s1 ="";
  String s2 = nr_pools;
  String s3 ="";
  int tseqlength =tseq_.length();
  String s4 =tseqlength;
  
  if (Ac_arma.n_rows != Ac_arma.n_cols)  s1 +=(" Matrix A must be quadratic! ");
  if (Ag_arma.n_rows != Ag_arma.n_cols)  s1 +=(" Matrix A must be quadratic! ");
  if (C0_arma.n_elem != nr_pools) s1 +=(" Length of C0 must match number of pools! ");   
  if (Cin_arma.n_cols != nr_pools) s1 +=(" Number of Columns of C input data frame must match number of pools! ");
  if (xi_arma.n_cols != nr_pools) s1 +=(" Number of Columns of Xi input data frame must match number of pools! ");    
  if (strlen(s1.get_cstring())>10)
  {
    s1 +=("Number of pools defined as: ");
    s1 +=s2;   
  }    
  if (Cin_arma.n_rows != tseq_.length()) s3 +=(" Number of Rows of C input data frame must match number of time steps! ");
  if (xi_arma.n_rows != tseq_.length()) s3 +=(" Number of Rows of Xi input data frame must match number of time steps! ");  
  if (strlen(s3.get_cstring())>10)
  {
    s3 +=("Number of time steps defined as: ");
    s3 +=s4;   
  }          
  s1+=s3;
  if (strlen(s1.get_cstring())>10) throw exception(s1.get_cstring());
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //main part
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int nsteps=tseq_.length();
  NumericMatrix cpools(nsteps,nr_pools); 
  arma::mat cpools_arma = as<arma::mat>(cpools);
  
  arma::mat cpools_sub = arma::zeros(1,nr_pools); 
  
  arma::mat K1_sub(1,nr_pools); 
  arma::mat K2_sub(1,nr_pools); 
  arma::mat K3_sub(1,nr_pools); 
  
  // start values for C
  cpools_arma.row(0)=arma::trans(C0_arma);
  
  // time loop
  for(int ts=1; ts<nsteps; ts++) 
  {  
    arma::rowvec xi_new_0 = xi_arma.row(ts-1);
    arma::rowvec xi_new_1 = xi_arma.row(ts);
    arma::rowvec K1(nr_pools);
    arma::rowvec K2(nr_pools);
    arma::rowvec K3(nr_pools);
    arma::rowvec K4(nr_pools);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //C-Pools calculation
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    arma::rowvec Cin_here=Cin_arma.row(ts-1);
    arma::rowvec Cin_1=Cin_here;
    arma::rowvec Cin_4=Cin_here;
    
    if(meanCin==TRUE) // alternative version, mathematically more consequent, meanCin!=TRUE makes sense for single input events
    {
      Cin_here = (Cin_here + Cin_arma.row(ts))/2;
      Cin_4=Cin_arma.row(ts);
    }
    
    bool negatives=TRUE;
    int substeps=1;
    
    while (negatives)
    {
      arma::mat add_pools_sub= arma::zeros(substeps,nr_pools); 
      cpools_sub=join_vert(cpools_arma.row(ts-1),add_pools_sub);
      bool anynegatives=false;
      arma::rowvec xi_new_0_ss = xi_new_0;
      
      for(int ss=1; ss<substeps+1; ss++)
      {
        if (ss>1) xi_new_0_ss=xi_new_1;
        if (( !is_crop_(ts) ) & (cpools_arma(ts-1,2)< aggr_sat)) 
        {
          A_arma = Ag_arma;
        }else{
          A_arma = Ac_arma;
        }
        K1=Cin_1/substeps + arma::trans((A_arma * arma::diagmat(xi_new_0_ss/timediv/substeps)  * arma::trans(cpools_sub.row(ss-1)))) ;   
        K2=Cin_here/substeps + arma::trans((A_arma * arma::diagmat((xi_new_0_ss+xi_new_1)/2/timediv/substeps)  * arma::trans(cpools_sub.row(ss-1)+0.5*K1))) ; 
        K3=Cin_here/substeps + arma::trans((A_arma * arma::diagmat((xi_new_0_ss+xi_new_1)/2/timediv/substeps)  * arma::trans(cpools_sub.row(ss-1)+0.5*K2))) ;  
        K4=Cin_4/substeps + arma::trans((A_arma * arma::diagmat(xi_new_1/timediv/substeps)  * arma::trans(cpools_sub.row(ss-1)+K3))) ;
        
        if (ss==1) 
        {
          K1_sub=K1;
          K2_sub=K2;
          K3_sub=K3;
        }
        
        if (ss>1)
        {
          K1_sub=join_vert(K1_sub, K1);
          K2_sub=join_vert(K2_sub, K2);
          K3_sub=join_vert(K3_sub, K3);
        }
        
        cpools_sub.row(ss)=(K1+ 2*K2 + 2*K3 +K4)/6+cpools_sub.row(ss-1);
        
        if (any(cpools_sub.row(ss)<0))
        {
          anynegatives=TRUE;
          break;
        }
        
      }
      
      if (anynegatives==TRUE)
      {
        substeps=substeps*2;
        Rcpp::Rcout << " substeps: "<< substeps << std::endl;
      }
      if (anynegatives==FALSE) negatives=false;
    }
    
    cpools_arma.row(ts)=cpools_sub.row(substeps);
  }
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //create output
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //define list output names
  CharacterVector list_names = {"C"};
  
  CharacterVector pool_i_names;
  
  for(int i=0; i<nr_pools; i++) 
  { 
    String pool_i_name("pool_");
    pool_i_name+=i+1;
    pool_i_names.push_back(pool_i_name);
  }
  
  //define length of output list
  int listlength = 1; //c always out
 
  List listout (listlength); 
  
  cpools=wrap(cpools_arma);
  colnames(cpools)=pool_i_names;
  
  listout(0)=cpools;
  
  
  
  listout.names()=list_names;
  //  readline("");
  //if (listlength<2) Rcout << "Ready. Produced list of length " << listlength << ". \n";
  return listout;
}

