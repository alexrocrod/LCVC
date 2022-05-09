% LCVC
% Módulo 3
% 16 de Maio de 2022
% Alexandre Rodrigues, 92993
% Anastasiya Lykholat, 93436
% Hugo Martins, 93247

% function simulation_SIR

clear all
close all
clc
format short

global alpha betaA betaP theta1 theta2 theta3 omega gamma epsilon miu miuA miuH sigma zeta niu; 

% - parameters --
alpha=0.2;
betaA=0.000273; 
betaP=0.000777; 
theta1=0.0003;
epsilon=1.5;
miu=0.00868; 
miuA=0.00775;   
miuH=0.0271;
gamma=0.00744;   
theta2=3*theta1; 
sigma=0.7;
zeta=0.0214;
theta3=16*theta1; 
niu=0.0155;
omega=0.0000000001;


% - final time -- 
T = 100; 

% -- initial conditions --
% P0=0.0553;
% A0=0.00148;
% H0=0.0003;
% R0=0.0097;
% S0=1-P0-A0-H0-R0;

% check R0
P0=0.12;
A0=0.0000000000000001;
H0=0.0000000000000001;
R0=0;
S0=1-P0-A0-H0-R0;

x0 = [S0; P0; A0; H0; R0]; 


%% -- solve system of ODE's ----

options = odeset('AbsTol',1e-12,'RelTol',1e-12) ;
[t, z] = ode45(@sys, [0 T], x0, options);


% -------- figures ----------------
%% 
figure; 
plot(t, z(:, 1),'-b', t, z(:, 2),'--r', t, z(:, 3), '.-g',t, z(:, 4),'-*y', t, z(:, 5),'-|k', 'LineWidth',1.5); 
xlabel('time');  ylabel('Number individuals'); 
legend('S', 'P', 'A', 'H', 'R')

%% 
figure; 
plot(t, z(:, 2),'--r', t, z(:, 3), '.-g',t, z(:, 4),'-*y', t, z(:, 5),'-|k', 'LineWidth',1.5); 
xlabel('time');  ylabel('Number individuals'); 
legend('P', 'A', 'H', 'R')

%%
figure; 
plot(t, z(:, 3), '.-g',t, z(:, 4),'-*y', t, z(:, 5),'-|k', 'LineWidth',1.5); 
xlabel('time');  ylabel('Number individuals'); 
legend('A', 'H', 'R')

%% Making sure S+P+A+H+R=1

erro1 = mean(sum(z,2) - ones(size(z,1),1))
  
%%

% -- System SIR -------------------

function zdot=sys(t,z)

global alpha betaA betaP theta1 theta2 theta3 omega gamma epsilon miu miuA miuH sigma zeta niu; 

S=z(1); P=z(2); A=z(3); H = z(4); R = z(5);

zdot = [ - alpha.* S - betaA.* S.* A - betaP.* S.* P - theta1.* S.* H + epsilon.* P + miu.*(P+A+H+R) + miuA.* A + miuH.* H;
         alpha.* S - epsilon.* P - gamma.* P - theta2.* P.* H - miu.* P;
         gamma.* P + sigma.* R.* A./(A+H+omega) + betaA.* S.* A + betaP.* S.* P - zeta.* A - theta3.* A.* H - (miu+miuA).* A;
         theta1.* S.* H + theta2.* P.* H + theta3.* A.* H + sigma.* R.* H./(A+H+omega) - niu.* H - (miu+miuH).* H;
         zeta.* A + niu.* H - sigma.* R.* A./(A+H+omega) - sigma.* R.* H./(A+H+omega) - miu.* R] ; 

end
%% controlo (u)  substituir o alpha (passsagem de suscetiveis para prescritos)


% zdot = [ - alpha.* S - betaA.* S.* A - betaP.* S.* P - theta1.* S.* H + epsilon.* P + miu.*(P+A+H+R) + miuA.* A + miuH.* H;
%          alpha.* S - epsilon.* P - gamma.* P - theta2.* P.* H - miu.* P;
%          gamma.* P + sigma.* R.* A./(A+H+omega) + betaA.* S.* A + betaP.* S.* P - zeta.* A - theta3.* A.* H - (miu+miuA).* A;
%          theta1.* S.* H + theta2.* P.* H + theta3.* A.* H + sigma.* R.* H./(A+H+omega) - niu.* H - (miu+miuH).* H;
%          zeta.* A + niu.* H - sigma.* R.* A./(A+H+omega) - sigma.* R.* H./(A+H+omega) - miu.* R] ; 


% maximimar prescricoes P minimizando o vicio A e H

%% controlo (u) (adicionar ao niu) (passsagem de Heroina (H) para recuperados (R)) novo tratamento


% zdot = [ - alpha.* S - betaA.* S.* A - betaP.* S.* P - theta1.* S.* H + epsilon.* P + miu.*(P+A+H+R) + miuA.* A + miuH.* H;
%          alpha.* S - epsilon.* P - gamma.* P - theta2.* P.* H - miu.* P;
%          gamma.* P + sigma.* R.* A./(A+H+omega) + betaA.* S.* A + betaP.* S.* P - zeta.* A - theta3.* A.* H - (miu+miuA).* A;
%          theta1.* S.* H + theta2.* P.* H + theta3.* A.* H + sigma.* R.* H./(A+H+omega) - niu.* H - (miu+miuH).* H;
%          zeta.* A + niu.* H - sigma.* R.* A./(A+H+omega) - sigma.* R.* H./(A+H+omega) - miu.* R] ; 

% zdot = [ - alpha.* S - betaA.* S.* A - betaP.* S.* P - theta1.* S.* H + epsilon.* P + miu.*(P+A+H+R) + miuA.* A + miuH.* H;
%          alpha.* S - epsilon.* P - gamma.* P - theta2.* P.* H - miu.* P;
%          gamma.* P + sigma.* R.* A./(A+H+omega) + betaA.* S.* A + betaP.* S.* P - zeta.* A - theta3.* A.* H - (miu+miuA).* A;
%          theta1.* S.* H + theta2.* P.* H + theta3.* A.* H + sigma.* R.* H./(A+H+omega) - niu.* H - (miu+miuH).* H -u*H; 
%          zeta.* A + niu.* H - sigma.* R.* A./(A+H+omega) - sigma.* R.* H./(A+H+omega) - miu.* R + u*H] ; 


% maximimar recuperados R minimizando tratamentos u (reduz custo)


