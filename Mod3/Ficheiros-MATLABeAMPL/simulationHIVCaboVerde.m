
function simulation  

% model HIV/AIDS - Cabo Verde 

clear all ; close all ; clc ; format short ;

global Lambda mu beta omega alpha rho phi d etaC etaA; 

% -  initial conditions ---

S0 = 338984-61; I0 = 61; C0 = 0; A0 = 0; 
N0 = S0 + I0 + C0 + A0; 

x0 = [S0; I0; C0; A0]; 

% -- final time in years --------
tf = 27;


% -- Parameters ------

mu = 1/69.54; 
Lambda = 2.2*1/69.54*N0; % =  1.0724e+004
beta = 0.752; 
etaC = 0.015;
etaA = 1.3;
omega = 0.09; 
rho = 0.1; 
phi = 1;
alpha = 0.33; 
d = 1;   

% ----------- basic reproduction number ---------------------
C1 = alpha + mu + d;
C2 = omega + mu; 
C3 = rho + phi + mu;

NumerR0 = beta*(C2*(C1 + etaA*rho) + etaC*phi*C1); 

DenomR0 = (mu*(C2*(rho+ C1)+C1*phi + d*rho) + rho*omega*d); 


R0 = NumerR0/DenomR0
     

% -------------------- solve ODE ---------------------
options = odeset('AbsTol',1e-12,'RelTol',1e-12) ;
[t, z] = ode45(@sysHIV,[0 tf],[x0], options);

S = z(:, 1); 
I = z(:, 2);
C = z(:, 3);
A = z(:, 4);

% -- real data from Cabo Verde 1986-2014  -----

totalHIVAIDS=[61,107,16,211,244,303,337,358,395,432,471,560,660,779,913,...
1064,1233,1493,1716,2015,2334,2610,2929,3340,3739,4090,4537,4946]; 

timeyears=0:1:27;


% --------------------- Figures --------------------------------------

% ---- compare real data with HIV/AIDS model 

figure
hold on
plot(timeyears,totalHIVAIDS,'red*','LineWidth',2.5)
plot(t, I + C + A + mu*(I + C) + (mu+d)*A,'b', 'LineWidth',1.5);
xlabel('time (years)'); ylabel('cumulative HIV and AIDS cases'); 
axis([0 28 0 5500]); 
legend('real data', 'model'); 
title('\eta_A = 1.3 and \eta_C = 0.015'); 
hold off


%--------------------------


function zdot=sysHIV(t,z)

global Lambda etaC etaA mu beta omega alpha rho phi d;

% x = [S=x1; I=x2; C=x3; A=x4]; 

x1=z(1); x2=z(2); x3=z(3); x4=z(4); 

zdot = [ Lambda - beta/(x1+x2+x3+x4).*(x2 + etaC.*x3 + etaA.*x4).*x1 - mu.*x1
         beta/(x1+x2+x3+x4).*(x2 + etaC.*x3 + etaA.*x4).*x1 - (rho + phi + mu).*x2 + alpha.*x4 + omega.*x3  
         phi.*x2 - (omega + mu).*x3
         rho.*x2 - (alpha + mu + d).*x4] ;  % system of ode's



