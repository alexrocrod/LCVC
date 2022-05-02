

####################################################
#                                                  #
#        Optimal control - covidPT - SAIRP         #
#                                                  #
#                                                  #
####################################################


 param tf   := 120;
 param n    := 1200 ;
 param h    := tf/n;
 
 param u1max := 0.5 ; # podemos considerar outros valores 0 < umax < 1
 
 param K1 = 100; 
 param K2 = 1; 
 
 ## initial values of state variables
 
 param N0 := 10295909; 
 
 ## cada classe representa a fracao em relacao a populacao total
 
 param S0 := 10295894/N0;
 param A0 := (2/0.15)/N0;
 param I0 := 2/N0 ;
 param P0 := 0/N0 ;
 
 
##  parameters 

 param beta := 1.464;
 param theta := 1;
 param p := 0.675;
 param phi := 1/12;
 
 param v := 1/1;
 param q := 0.15;
 param nu := v*q; 
 
 param delta := 1/27; 
 
 param m := 0.09; 
 param w := 1/45;
 
                        
##  state variables

 var int {i in 0..n};
 s.t. iv_int : int[0] = 0 ;

 var S {i in 0..n};
 s.t. iv_S : S[0] = S0 ;
 s.t. mu_V {i in 0..n} : S[i] >= 0 ;
 
 var A {i in 0..n};
 s.t. iv_A : A[0] = A0 ;
 s.t. mu_A {i in 0..n} : A[i] >= 0 ;
 
 var I {i in 0..n};
 s.t. iv_I : I[0] = I0 ;
 s.t. mu_I {i in 0..n} : I[i] >= 0 ;
 s.t. sc_I {i in 0..n} : I[i] <= 2/3*2.5*10^(-3) ;
 
 var P {i in 0..n};  
 s.t. iv_P : P[0] = P0 ;
 s.t. mu_P {i in 0..n} : P[i] >= 0 ;
  

 
## control variables

 var u1 {i in 0..n}   ;
 s.t. mu_c1 {i in 0..n} : 0 <= u1[i] <= u1max  ;


## right hand sides of ODEs

 var fint {i in 0..n} =  K1*I[i] - K2*u1[i];

 var fS {i in 0..n} = - beta*(theta*A[i]+I[i])*(1-p)*S[i] - phi*p*S[i] + w*u1[i]*P[i];
 
 var fA {i in 0..n} = beta*(theta*A[i]+I[i])*(1-p)*S[i] - nu*A[i];
 
 var fI {i in 0..n} = nu*A[i] - delta*I[i] ;
   
 var fP {i in 0..n} = phi*p*S[i] - w*u1[i]*P[i] ;
 
                                       

##  objective functional

 minimize OBJ :  int[n] ;

## Trapezoidal method for ODEs


 s.t. l_int {i in 0..n-1} : int[i+1] = int[i] + 0.5*h * (fint[i]+fint[i+1]) ;

 s.t. l_S {i in 0..n-1}  : S[i+1] = S[i] + 0.5*h * (fS[i]+fS[i+1])  ;
 s.t. l_A {i in 0..n-1}  : A[i+1] = A[i] + 0.5*h * (fA[i]+fA[i+1])  ; 
 s.t. l_I {i in 0..n-1}  : I[i+1] = I[i] + 0.5*h * (fI[i]+fI[i+1])  ;
 s.t. l_P {i in 0..n-1}  : P[i+1] = P[i] + 0.5*h * (fP[i]+fP[i+1])  ;


#####  NEOS SOLVER IPOPT   ##########
## - https://neos-server.org/neos/solvers/nco:Ipopt/AMPL.html

option abs_boundtol 1;

option solver ipopt;
option ipopt_options " max_iter=2000 acceptable_tol=1e-10 "; 

solve;


###### OUTPUT SCREEN

display OBJ; 

printf{i in 0..n}: "%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
i*tf/n, S[i], A[i], I[i], P[i], u1[i];

## end; 




