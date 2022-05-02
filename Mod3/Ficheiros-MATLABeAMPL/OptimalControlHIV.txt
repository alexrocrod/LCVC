
####################################################
#                                                  #
#        Optimal control of HIV transmission       #
#           Normalized  SICA model                 #
#                                                  #
####################################################

 param tf   := 20;
 param n    := 2000;
 param h    := tf/n;
 
 param umax = 0.5 ;

 
 ## initial values of state variables
 
 param S0 := 0.6 ;
 param I0 := 0.2 ;
 param C0 := 0.1 ;
 param A0 := 0.1 ;
 

##  parameters 

 param mu := 1/69.54; 
 param b := 2.1*mu; 
 param beta := 1.9;
 param omega := 0.09; 
 param rho := 0.1; 
 param phi := 1;
 param alpha := 0.33; 
 param d := 1;   
 param etaC := 0.015; 
 param etaA := 1.3; 


##  state variables

 var int {i in 0..n};
 s.t. iv_int : int[0] = 0 ;

 var S {i in 0..n};
 s.t. iv_s : S[0] = S0 ;
 
 var I {i in 0..n};
 s.t. iv_I : I[0] = I0 ;
 
 var C {i in 0..n};  
 s.t. iv_C : C[0] = C0 ;
 
 var A {i in 0..n};  
 s.t. iv_A : A[0] = A0 ;
 
 
## control variables

 var u {i in 0..n} := umax  ;
 s.t. mu_c {i in 0..n} : 0 <= u[i] <= umax  ;



## right hand sides of ODEs

 var fint {i in 0..n} = S[i] - I[i] - u[i]*u[i];

 var fS {i in 0..n} = b - b*S[i] - (1 - u[i])*beta*(I[i] + etaC*C[i] + etaA*A[i])*S[i] + d*A[i]*S[i];
 
 var fI {i in 0..n} = (1-u[i])*beta*(I[i] + etaC*C[i] + etaA*A[i])*S[i] - (rho + phi + b)*I[i] + alpha*A[i]  + omega*C[i] + d*A[i]*I[i]  ;

 var fC {i in 0..n} = phi*I[i] - (omega + b)*C[i] +  d*A[i]*C[i] ; 
   
 var fA {i in 0..n} = rho*I[i] - (alpha + b + d)*A[i] + d*A[i]*A[i] ;
                                         

##  objective functional

 maximize OBJ :  int[n] ;

## EULER method for ODEs

 s.t. l_int {i in 0..n-1} : int[i+1] = int[i] + h * fint[i] ;

 
 ## Implicit EULER method for ODEs

 s.t. l_S {i in 0..n-1}  : S[i+1] = S[i] + 0.5*h * (fS[i]+fS[i+1])  ;
 s.t. l_I {i in 0..n-1}  : I[i+1] = I[i] + 0.5*h * (fI[i]+fI[i+1])   ;
 s.t. l_C {i in 0..n-1}  : C[i+1] = C[i] + 0.5*h * (fC[i]+fC[i+1])  ;
 s.t. l_A {i in 0..n-1}  : A[i+1] = A[i] + 0.5*h * (fA[i]+fA[i+1])  ;
     

#####  NEOS SOLVER IPOPT   ##########
## - https://neos-server.org/neos/solvers/nco:Ipopt/AMPL.html

option abs_boundtol 1;

option solver ipopt;
option ipopt_options " max_iter=2000 acceptable_tol=1e-10 "; 

solve;


###### OUTPUT SCREEN

display OBJ; 

printf{i in 0..n}: "%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
i*tf/n, S[i], I[i], C[i], A[i], u[i];

end; 



