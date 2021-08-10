
####################################################
#                                                  #
#        Optimal control of a case study           #
#           of Normalized SEIR model               #
#                                                  #
####################################################

 param tf   := 2050-2014;
 param n    := 500;
 param h    := tf/n;
 
 param umax = 0.9 ;
 
 ## initial values of state variables
 
 param S0 := 0.88 ;
 param E0 := 0.07 ;
 param I0 := 0.05 ;
 param R0 := 0 ;
 

##  parameters 

 param mu := 0.1;  
 param beta := 1.9;
 param gamma := 0.1887; 
 param A := 0.3;


##  state variables

 var J {i in 0..n};
 s.t. iv_J : J[0] = 0 ;

 var S {i in 0..n};
 s.t. iv_s : S[0] = S0 ;
 
 var E {i in 0..n};  
 s.t. iv_E : E[0] = E0 ;

 var I {i in 0..n};
 s.t. iv_I : I[0] = I0 ;
 
 var R {i in 0..n};  
 s.t. iv_R : R[0] = R0 ;
 
 
## control variables

 var u {i in 0..n} := umax  ;
 s.t. mu_c {i in 0..n} : 0 <= u[i] <= umax  ;

## right hand sides of ODEs

 var fJ {i in 0..n} = I[i] + u[i]*u[i]; 

 var fS {i in 0..n} = - beta*S[i]*I[i] - u[i]*S[i];
 
 var fE {i in 0..n} = beta*S[i]*I[i] - gamma*E[i]  ;

 var fI {i in 0..n} = gamma*E[i] - mu*I[i] ; 
   
 var fR {i in 0..n} = mu*I[i] + u[i]*S[i];
                                         

##  objective functional

 minimize OBJ :  J[n] ;


## EULER method for ODEs

 s.t. l_J {i in 0..n-1} : J[i+1] = J[i] + h * fJ[i] ;
 

## Implicit EULER method for ODEs

 s.t. l_S {i in 0..n-1}  : S[i+1] = S[i] + 0.5*h * (fS[i]+fS[i+1])  ;
 s.t. l_E {i in 0..n-1}  : E[i+1] = E[i] + 0.5*h * (fE[i]+fE[i+1])  ;
 s.t. l_I {i in 0..n-1}  : I[i+1] = I[i] + 0.5*h * (fI[i]+fI[i+1])   ;
 s.t. l_R {i in 0..n-1}  : R[i+1] = R[i] + 0.5*h * (fR[i]+fR[i+1])  ;
     

#####  NEOS SOLVER IPOPT   ##########
## - https://neos-server.org/neos/solvers/nco:Ipopt/AMPL.html

option abs_boundtol 1;

option solver ipopt;
option ipopt_options " max_iter=500 acceptable_tol=1e-10 "; 

solve;


###### OUTPUT SCREEN

display OBJ; 

printf{i in 0..n}: "%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
i*tf/n, S[i],  E[i], I[i], R[i], u[i];

end; 




