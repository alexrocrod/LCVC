function simulation_SIR

clear all ; close all ; clc ; format short ;

global beta gamma; 

% - parameters --
beta = 0.0055; % 0.0055
gamma = 0.33;  % 0.33


% - final time -- 
T = 100; 

% -- initial conditions --

S0 = 290; 
I0 = 10; 
R0 = 0; 

x0 = [S0; I0; R0]; 

% -- solve system of ODE's ----

options = odeset('AbsTol',1e-12,'RelTol',1e-12) ;
[t, z] = ode45(@sys, [0 T], x0, options);


% -------- figures ----------------

figure; 
plot(t, z(:, 1),'-b', t, z(:, 2),'--r', t, z(:, 3), '.-g', 'LineWidth',1.5); 
xlabel('time');  ylabel('Number individuals'); 
legend('S', 'I', 'R')



% -- System SIR -------------------

function zdot=sys(t,z)

global beta gamma;

% x = [S=x1; I=x2; R=x3]; 

x1=z(1); x2=z(2); x3=z(3); 

zdot = [ - beta.*x1.*x2 
         beta.*x1.*x2 - gamma.*x2  
         gamma.*x2] ; 
