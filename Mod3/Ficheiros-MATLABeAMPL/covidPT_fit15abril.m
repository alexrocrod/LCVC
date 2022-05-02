function covidPT_fit15abril

clear all; close all; clc; format long;

global beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 beta9 theta sigma p1 p2 p3 p4 p5 p6 p7 p8 p9 phi w nu delta m1 m2 m3 m4 m5 m6 m7 m8 m9 Lambda mu; 

%-----------------------------------------------
%           Initial conditions  
%------------------------------------------------
%  
N0 = 10295894 + 2 + 2/0.15;

S0 = 10295894; 
I0 = 2; 
A0 = (2/0.15);
R0 = 0; 
P0 = 0; 

%---------------------------------------
%               real data
%---------------------------------------

M=xlsread('data-pt-covid19-15abril');

tf = length(M) %
time=1:1:tf;   % consideramos 410 dias

% ---- consider n subintervals of time --- 

tf1 = 73;  % 2 marco, 2020 - 13 maio
tf2 = 90;  %  30 maio 
tf3 = 130; %  9 julho
tf4 = 163; %  11 agosto
tf5 = 200; %  17 setembro
tf6 = 253; %  9 novembro
tf7 = 304; %  30 dezembro
tf8 = 329; %  24 janeiro
% tf = 15 abril

ativos = M(1:tf,2)

% % ---- figure active cases 
% figure; 
% hold on
% plot(time, ativos, 'r*'); 
% x=xlabel('Time (days)');
% y=ylabel('Active infected individuals - COVID-19');
% t=title('COVID-19 Portuguese data: 02-03-202 -- 14-05-2021', 'FontSize', 12); 
% set(y, 'FontName', 'Arial', 'FontSize', 12);
% set(x, 'FontName', 'Arial', 'FontSize', 12);
% set(t, 'FontName', 'Arial', 'FontSize', 12);


%  ---------------- Parameters - fixed ------------------------
mu = 1/(81*365); 
Lambda = (0.19/100*N0)/365;
theta = 1;
sigma = 0;
                                    
phi = 1/12;                                                               
v = 1/1;                                                                   
q = 0.15;                                                                 
nu = v*q;
delta = 1/27;                                                     
w =  1/45;  

% -- beta para os n intervalos de tempo 
 
beta1 = 1.502; 
beta2 = 0.6; 
beta3 = 1.24; 
beta4 = 0.936;
beta5 = 1.531;
beta6 = 0.886;
beta7 = 0.250;
beta8 = 0.793;
beta9 = 0.1;  

% -- p para os n intervalos de tempo - https://www.pse.pt/evolucao-confinamento-mobilidade/
p1 =  0.675; %(4950354+2000000)/N0;     2 March, 2020 - 13 May, 2020                                            % p: percentagem de pessoas que se movem de S para Q.
p2 = 0.650;  % 13 May - 30 May
p3 =  0.58;  %0.580; % 30 May - 9 July
p4 =  0.610; % 9 July - 11 August 
p5 = 0.580 ; % 11 August - 17 September
p6 = 0.290 ; % 17 September - 9 November
p7 = 0.370 ; % 9 November - 30 December
p8 = 0.370 ; % 30 December, 2020 - 24 January, 2021
p9 = 0.550 ; % 24 January - 15 April 


% -- m para os n intervalos de tempo
m1 = 0.066; 
m2 = 0.090;
m3 = 0.180;  
m4 = 0.160; 
m5 = 0.170; 
m6 = 0.140; 
m7 = 0.379; 
m8 = 0.090; 
m9 =  0.09; 


% ---- resolver ode entre 2 de marco 2020 e 15 abril 2021
options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[tn1, xn1] = ode45(@sys_norm1,[0 tf1],[S0; A0; I0; R0; P0],options);
% ---- resolver ode entre 13 maio e 30 maio 
[tn2, xn2] = ode45(@sys_norm2,[tf1 tf2],[xn1(end,1); xn1(end,2); xn1(end,3); xn1(end,4); xn1(end,5)],options);
% --- resolver ode entre 30 maio a 9 julho 
[tn3, xn3] = ode45(@sys_norm3,[tf2 tf3],[xn2(end,1); xn2(end,2); xn2(end,3); xn2(end,4); xn2(end,5)],options);
% ---- resolver ode entre 9 julho a 11 de agosto 
[tn4, xn4] = ode45(@sys_norm4,[tf3 tf4],[xn3(end,1); xn3(end,2); xn3(end,3); xn3(end,4); xn3(end,5)],options);
% ---- resolver ode 11 de agosto e 17 setembro 
[tn5, xn5] = ode45(@sys_norm5,[tf4 tf5],[xn4(end,1); xn4(end,2); xn4(end,3); xn4(end,4); xn4(end,5)],options);
% ---- resolver ode 17 setembro e 9 novembro 
[tn6, xn6] = ode45(@sys_norm6,[tf5 tf6],[xn5(end,1); xn5(end,2); xn5(end,3); xn5(end,4); xn5(end,5)],options);
% ---- resolver ode 9 novembro e 30 dezembro, 2020
[tn7, xn7] = ode45(@sys_norm7,[tf6 tf7],[xn6(end,1); xn6(end,2); xn6(end,3); xn6(end,4); xn6(end,5)],options);
% ---- resolver ode 30 dezembro2020 e 24 janeiro 2021 
[tn8, xn8] = ode45(@sys_norm8,[tf7 tf8],[xn7(end,1); xn7(end,2); xn7(end,3); xn7(end,4); xn7(end,5)],options);
% ---- resolver ode 24 janeiro ate 15 abril = tf
[tn9, xn9] = ode45(@sys_norm9,[tf8 tf],[xn8(end,1); xn8(end,2); xn8(end,3); xn8(end,4); xn8(end,5)],options);

% ---------- Figure -------------

figure; 
hold on;
plot(tn1, xn1(:,3),'-r','LineWidth',3.5);
plot(tn2, xn2(:,3),'-b','LineWidth',3.5);
plot(tn3, xn3(:,3),'-c','LineWidth',3.5);
plot(tn4, xn4(:,3),'-m','LineWidth',3.5);
plot(tn5, xn5(:,3),'-y','LineWidth',3.5);
plot(tn6, xn6(:,3),'Color',[0.9290 .6940 0.1250],'LineStyle','-','LineWidth',3.5);
plot(tn7, xn7(:,3),'Color',[0.4660 .6740 0.1880],'LineStyle','-','LineWidth',3.5);
plot(tn8, xn8(:,3),'Color',[0.6350 .0780 0.1840],'LineStyle','-','LineWidth',3.5);
plot(tn9, xn9(:,3),'-g','LineWidth',3.5);
plot(time, M(time,2), ':k','LineWidth',2); 
l=legend('March 2, 2020 - May 13', 'May 13 - May 30, 2020', 'May 30 - July 9', 'July 9 - August 11', 'August 11 - September 17',...
    'September 17 - November 9', 'November 9 - December 30', 'December 30, 2020 - January 24, 2021', 'January 24 - April 15','Real data', 'Location', 'southeastoutside', 'boxoff');
y1 = ylim; % current y-axis limits
plot([tf1 tf1],[y1(1) y1(2)],'r:','LineWidth',1.5);
plot([tf2 tf2],[y1(1) y1(2)], 'b:','LineWidth',1.5);
plot([tf3 tf3],[y1(1) y1(2)],'c:','LineWidth',1.5);
plot([tf4 tf4],[y1(1) y1(2)],'m:','LineWidth',1.5);
plot([tf5 tf5],[y1(1) y1(2)], 'y:','LineWidth',1.5);
plot([tf6 tf6],[y1(1) y1(2)], 'Color',[0.9290 .6940 0.1250],'LineStyle',':','LineWidth',1.5);
plot([tf7 tf7],[y1(1) y1(2)], 'Color',[0.4660 .6740 0.1880],'LineStyle',':','LineWidth',1.5);
plot([tf8 tf8],[y1(1) y1(2)], 'Color',[0.6350 .0780 0.1840],'LineStyle',':','LineWidth',1.5);
x=xlabel('Time (days)');
y=ylabel('Active infected - I')
t=title('Fit COVID-19 Portuguese data', 'FontSize', 12); 
set(t, 'FontName', 'Arial', 'FontSize', 12);
set(l, 'FontName', 'Arial', 'FontSize', 12);
set(y, 'FontName', 'Arial', 'FontSize', 12);
set(x, 'FontName', 'Arial', 'FontSize', 12);


% --------------- 

function xdot1 = sys_norm1(t,x)
global  beta1 theta p1 phi w nu delta m1 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];

x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
C=beta1*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot1 = [ Lambda - C.*(1-p1)*x1 - phi*p1*x1 + w*m1*x5 - mu*x1  % S
         C.*(1-p1)*x1 - nu*x2  - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                       % I 
         delta*x3 - mu*x4                               % R
         phi*p1*x1 - w*m1*x5 - mu*x5 ];                  % P


% --------------- 

function xdot2 = sys_norm2(t,x)
global  beta2 theta p2 phi w nu delta m2 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];
x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 

C=beta2*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot2 = [ Lambda - C.*(1-p2)*x1 - phi*p2*x1 + w*m2*x5 - mu*x1  % S
         C.*(1-p2)*x1 - nu*x2 - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                      % I 
         delta*x3 - mu*x4                              % R
         phi*p2*x1 - w*m2*x5 - mu*x5];                 % P

% --------------- 

function xdot3 = sys_norm3(t,x)
global  beta3 theta p3 phi w nu delta m3 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];

x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
C=beta3*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot3 = [ Lambda - C.*(1-p3)*x1 - phi*p3*x1 + w*m3*x5 - mu*x1  % S
         C.*(1-p3)*x1 - nu*x2 - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                      % I 
         delta*x3  - mu*x4                             % R
         phi*p3*x1 - w*m3*x5 - mu*x5 ];                 % P
     
% --------------- 

function xdot4 = sys_norm4(t,x)
global  beta4 theta p4 phi w nu delta m4 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];

x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
C=beta4*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot4 = [ Lambda - C.*(1-p4)*x1 - phi*p4*x1 + w*m4*x5 - mu*x1  % S
         C.*(1-p4)*x1 - nu*x2 - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                      % I 
         delta*x3  - mu*x4                             % R
         phi*p4*x1 - w*m4*x5 - mu*x5 ];                 % P
     
% --------------- 

function xdot5 = sys_norm5(t,x)
global  beta5 theta p5 phi w nu delta m5 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];

x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
C=beta5*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot5 = [ Lambda - C.*(1-p5)*x1 - phi*p5*x1 + w*m5*x5 - mu*x1  % S
         C.*(1-p5)*x1 - nu*x2 - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                      % I 
         delta*x3  - mu*x4                             % R
         phi*p5*x1 - w*m5*x5 - mu*x5 ];                 % P
     
% --------------- 

function xdot6 = sys_norm6(t,x)
global  beta6 theta p6 phi w nu delta m6 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];

x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
C=beta6*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot6 = [ Lambda - C.*(1-p6)*x1 - phi*p6*x1 + w*m6*x5 - mu*x1  % S
         C.*(1-p6)*x1 - nu*x2 - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                      % I 
         delta*x3  - mu*x4                             % R
         phi*p6*x1 - w*m6*x5 - mu*x5 ];                 % P
     
% --------------- 

function xdot7 = sys_norm7(t,x)
global  beta7 theta p7 phi w nu delta m7 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];

x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
C=beta7*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot7 = [ Lambda - C.*(1-p7)*x1 - phi*p7*x1 + w*m7*x5 - mu*x1  % S
         C.*(1-p7)*x1 - nu*x2 - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                      % I 
         delta*x3  - mu*x4                             % R
         phi*p7*x1 - w*m7*x5 - mu*x5 ];                 % P
     
% --------------- 

function xdot8 = sys_norm8(t,x)
global  beta8 theta p8 phi w nu delta m8 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];

x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
C=beta8*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot8 = [ Lambda - C.*(1-p8)*x1 - phi*p8*x1 + w*m8*x5 - mu*x1  % S
         C.*(1-p8)*x1 - nu*x2 - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                      % I 
         delta*x3  - mu*x4                             % R
         phi*p8*x1 - w*m8*x5 - mu*x5 ];                 % P
     
% --------------- 

function xdot9 = sys_norm9(t,x)
global  beta9 theta p9 phi w nu delta m9 Lambda mu;
% x = [S=x1; A=x2; I=x3; R=x4; P=x5];

x1=x(1); x2=x(2); x3=x(3); x4=x(4); x5=x(5); 
C=beta9*(theta*x2+x3)./(x1+x2+x3+x4+x5);

xdot9 = [ Lambda - C.*(1-p9)*x1 - phi*p9*x1 + w*m9*x5 - mu*x1  % S
         C.*(1-p9)*x1 - nu*x2 - mu*x2                  % A
         nu*x2 - delta*x3 - mu*x3                      % I 
         delta*x3  - mu*x4                             % R
         phi*p9*x1 - w*m9*x5 - mu*x5 ];                 % P
     

