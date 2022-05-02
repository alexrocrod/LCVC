
function simuHIVnormalized 

% -- HIV/AIDS - sica normalizado 

clear all ; close all ; clc ; format short ;

global b mu beta omega alpha rho phi d etaC etaA; 

% ----------------- Parameters -----------------------

mu = 1/69.54; 
b = 2.1*mu; 
beta = 1.6;
omega = 0.09; 
rho = 0.1; 
phi = 1;
alpha = 0.33; 
d = 1;   
etaC = 0.015; 
etaA = 1.3; 


% ------------- initial conditions ----------------------------------
s0 = 0.6; i0 = 0.2; c0 = 0.1; a0 = 0.1; 
x0 = [s0; i0; c0; a0]; 


% ------------------- final time ---------------------------
tf = 20; % (years) 

% -------------------- solve ODE ------------------------
options = odeset('AbsTol',1e-12,'RelTol',1e-12) ;
[t, z] = ode45(@sys,[0 tf],[x0], options);

s = z(:, 1); 
i = z(:, 2);
c = z(:, 3);
a = z(:, 4);

% --------------------- Figures --------------------------------------

figure
hold on
plot(t, s, t, i, t, c, t, a, 'LineWidth',1.5);
xlabel('time (years)'); ylabel('state variables'); 
axis([0 28 0 1]); 
legend('s', 'i', 'c', 'a'); 

figure; 
plot(t, s, t, i+c+a, 'LineWidth',1.5);
xlabel('time (years)'); ylabel('state variables'); 
axis([0 28 0 1]); 
legend('s', 'HIV-Infected + Chronic + AIDS'); 



% -------------------------------------------------

function zdot=sys(t,z)

global b mu beta omega alpha rho phi d etaC etaA;

% x = [s=x1; i=x2; c=x3; a=x4]; 

x1=z(1); x2=z(2); x3=z(3); x4=z(4); 

zdot = [ b  - b.*(x1) - beta.*(x2 + etaC.*x3 + etaA.*x4).*x1 + d.*x4.*x1
         beta.*(x2 + etaC.*x3 + etaA.*x4).*x1 - (rho + phi + b).*x2 + alpha.*x4  + omega.*x3 + d.*x4.*x2          
         phi.*x2 - (omega + b).*x3 +  d.*x4.*x3
         rho.*x2 - (alpha + b + d).*x4 + d.*x4.*x4] ;  % sistema de ode's
     
           
     
     
     
