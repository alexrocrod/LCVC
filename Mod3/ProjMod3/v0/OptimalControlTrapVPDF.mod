
####################################################
#                                                  #
#        Optimal control of HIV transmission       #
#           Normalized  SICA model                 #
#                                                  #
####################################################

 param tf   := 10;
 param n    := 1000;
 param h    := tf/n;
 
 param umax = 0.8;

 
 ## initial values of state variables
 
 param P0 := 0.095 ;
 param A0 := 0.0071;
 param H0 := 0.000465;
 param R0 := 0.0507 ;
 param S0 := 1-P0-A0-H0-R0 ;
 

##  parameters 

 param alpha := 0.27;
 param betaA := 0.000878; 
 param betaP := 0.0000654; 
 param theta1 := 0.222;
 param epsilon := 2.53;
 param miu := 0.0071;
 param miuA := 0.00883;
 param miuH := 0.0466;
 param gamma := 0.00505;
 param theta2 := 0.236;
 param sigma := 0.102;
 param zeta := 0.198;
 param theta3 := 19.7;
 param niu := 0.000531;
 param omega := 0.0000000001;


##  state variables

 var int {i in 0..n};
 s.t. iv_int : int[0] = 0 ;

 var S {i in 0..n};
 s.t. iv_s : S[0] = S0 ;

 var P {i in 0..n};
 s.t. iv_P : P[0] = P0 ;

 var A {i in 0..n};  
 s.t. iv_A : A[0] = A0 ;

 var H {i in 0..n};  
 s.t. iv_H : H[0] = H0 ;

 var R {i in 0..n};
 s.t. iv_R : R[0] = R0 ;


 
## control variables

 var u {i in 0..n} := umax  ;
 s.t. mu_c {i in 0..n} : 0 <= u[i] <= umax  ;



## right hand sides of ODEs

 var fint {i in 0..n} = R[i] - u[i]*u[i];

 var fS {i in 0..n} = - alpha*S[i] - betaA*S[i]*A[i] - betaP*S[i]*P[i] - theta1*S[i]*H[i] + epsilon*P[i] + miu*(P[i]+A[i]+H[i]+R[i]) + miuA*A[i] + miuH*H[i] ;
 
 var fP {i in 0..n} = alpha*S[i] - epsilon*P[i] - gamma*P[i] - theta2*P[i]*H[i] - miu*P[i];

 var fA {i in 0..n} = gamma*P[i] + sigma*R[i]*A[i]/(A[i]+H[i]+omega) + betaA*S[i]*A[i] + betaP*S[i]*P[i] - zeta*A[i] - theta3*A[i]*H[i] - (miu+miuA)*A[i];

 var fH {i in 0..n} = theta1*S[i]*H[i] + theta2*P[i]*H[i] + theta3*A[i]*H[i] + sigma*R[i]*H[i]/(A[i]+H[i]+omega) - niu*H[i] - (miu+miuH)*H[i] -u[i]*H[i];
 
 var fR {i in 0..n} = zeta*A[i] + niu*H[i] - sigma*R[i]*A[i]/(A[i]+H[i]+omega) - sigma*R[i]*H[i]/(A[i]+H[i]+omega) - miu*R[i] + u[i]*H[i] ;
 
                                         

##  objective functional

 maximize OBJ :  int[n] ;

## Trapezoidal method for ODEs

 s.t. l_int {i in 0..n-1} : int[i+1] = int[i] + 0.5*h * (fint[i]+fint[i+1]) ;

 s.t. l_S {i in 0..n-1}  : S[i+1] = S[i] + 0.5*h * (fS[i]+fS[i+1])  ;
 s.t. l_P {i in 0..n-1}  : P[i+1] = P[i] + 0.5*h * (fP[i]+fP[i+1])   ;
 s.t. l_A {i in 0..n-1}  : A[i+1] = A[i] + 0.5*h * (fA[i]+fA[i+1])  ;
 s.t. l_H {i in 0..n-1}  : H[i+1] = H[i] + 0.5*h * (fH[i]+fH[i+1])  ;
 s.t. l_R {i in 0..n-1}  : R[i+1] = R[i] + 0.5*h * (fR[i]+fR[i+1])  ;
     

#####  NEOS SOLVER IPOPT   ##########
## - https://neos-server.org/neos/solvers/nco:Ipopt/AMPL.html

option abs_boundtol 1;

option solver ipopt;
option ipopt_options " max_iter=2000 acceptable_tol=1e-10 "; 

solve;


###### OUTPUT SCREEN

display OBJ; 

printf{i in 0..n}: "%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
i*tf/n, S[i], P[i], A[i], H[i], R[i], u[i];

end; 



