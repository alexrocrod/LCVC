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

%% Do PDF ParamNorm
alpha=0.27;
betaA=0.000878; 
betaP=0.0000654; 
theta1=0.222;
epsilon=2.53;
miu=0.0071; 
miuA=0.00883;   
miuH=0.0466;
gamma=0.00505;   
theta2=0.236; 
sigma=0.102;
zeta=0.198;
theta3=19.7; 
niu=0.000531;
omega=0.0000000001; %%% onde esta no pdf

%% Parametros simplificados
% alpha=0.2;
% betaA=0.0001; 
% betaP=0.00001; 
% theta1=0.2;
% epsilon=2.5;
% miu=0; 
% miuA=0;   
% miuH=0;
% gamma=0.005;   
% theta2=0.2; 
% sigma=0.1;
% zeta=0.2;
% theta3=20; 
% niu=0.0005;
% omega=0.0000000001;

%% Parametros muito maus ParamMau

% alpha=0.3;
% betaA=0.01; 
% betaP=0.001; 
% theta1=0.5;
% epsilon=1;
% miu=0; 
% miuA=0;   
% miuH=0;
% gamma=0.05;   
% theta2=0.5; 
% sigma=0.2;
% zeta=0.01;
% theta3=20; 
% niu=0.00001;
% omega=0.0000000001;


%% Parametros muito bons 
% alpha=0.01;
% betaA=0; 
% betaP=0; 
% theta1=0.01;
% epsilon=5;
% miu=0; 
% miuA=0;   
% miuH=0;
% gamma=0.01;   
% theta2=0.01; 
% sigma=0.01;
% zeta=0.01;
% theta3=1; 
% niu=0.1;
% omega=0.0000000001;

%% Parametros bons ParamBom
% alpha=0.1;
% betaA=0.0001; 
% betaP=0.00001; 
% theta1=0.1;
% epsilon=3;
% miu=0; 
% miuA=0;   
% miuH=0;
% gamma=0.01;   
% theta2=0.1; 
% sigma=0.05;
% zeta=0.1;
% theta3=5; 
% niu=0.05;
% omega=0.0000000001;

%% 
% PDF % usado no VarParam
% P0=0.095;
% A0=0.0071;
% H0=0.000465;
% R0=0.0507;

% Ma situacao inicial
% P0=0.05;
% A0=0.05;
% H0=0.005;
% R0=0.0507;

% Inicial Normal
P0=0.1;
A0=0.001;
H0=0.0001;
R0=0.001;

fig1 = "VariaInicial/Matlab/GeralPDF4.jpg";
fig2 = "VariaInicial/Matlab/I2PDF4.jpg";
fig3 = "VariaInicial/Matlab/PDF4.jpg";

%% 
% - final time -- 
T = 10; %100

S0=1-P0-A0-H0-R0;

x0 = [S0; P0; A0; H0; R0]; 


%% -- solve system of ODE's ----

options = odeset('AbsTol',1e-12,'RelTol',1e-12) ;
[t, z] = ode45(@sys, [0 T], x0, options);


%% Making sure S+P+A+H+R=1
erro1 = mean(sum(z,2) - ones(size(z,1),1));
fprintf("Erro:%d\n",erro1)  

%% «
fh = figure();
fh.WindowState = 'maximized';
% plot(t, z(:, 1),'-b', t, z(:, 2),'--r', t, z(:, 3), '.-g',t, z(:, 4),'-*y', t, z(:, 5),'-|k', 'LineWidth',1.5);
plot(t, z(:, 1), t, z(:, 2), t, z(:, 3),t, z(:, 4), t, z(:, 5), 'LineWidth',1.5); 
xlabel('time');  ylabel('Number individuals'); 
legend('S', 'P', 'A', 'H', 'R')
saveas(gcf,fig1)

%% 
fh = figure();
fh.WindowState = 'maximized'; 
% plot(t, z(:, 2),'--r', t, z(:, 3), '.-g',t, z(:, 4),'-*y', t, z(:, 5),'-|k', 'LineWidth',1.5);
plot(t, z(:, 2), t, z(:, 3),t, z(:, 4), t, z(:, 5), 'LineWidth',1.5);
xlabel('time');  ylabel('Number individuals'); 
legend('P', 'A', 'H', 'R')
saveas(gca,fig2)

% %%
% figure; 
% % plot(t, z(:, 3), '.-g',t, z(:, 4),'-*y', t, z(:, 5),'-|k', 'LineWidth',1.5); 
% plot(t, z(:, 3),t, z(:, 4), t, z(:, 5), 'LineWidth',1.5);
% xlabel('time');  ylabel('Number individuals'); 
% legend('A', 'H', 'R')


%% Plots
fh = figure();
fh.WindowState = 'maximized'; 
cols = 3;
rows = 2;

s = z(:, 1); 
p = z(:, 2); 
a = z(:, 3); 
h = z(:, 4);
r = z(:, 5); 

subplot(rows,cols,1)
plot(t,s,'LineWidth',1.5); 
title("Suscetiveis (S)") 

subplot(rows,cols,2)
plot(t,p,'LineWidth',1.5); 
title("Prescritos (P)") 

subplot(rows,cols,3)
plot(t,a,'LineWidth',1.5); 
title("Viciados em Opioides (A)") 

subplot(rows,cols,4)
plot(t,h,'LineWidth',1.5); 
title("Viciados em Heroina (H)") 

subplot(rows,cols,5)
plot(t,r,'LineWidth',1.5); 
title("Recuperados (R)") 

% subplot(rows,cols,6)
% plot(t,u,'LineWidth',1.5); 
% title("Controlo (u)") 

saveas(gca,fig3)


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


