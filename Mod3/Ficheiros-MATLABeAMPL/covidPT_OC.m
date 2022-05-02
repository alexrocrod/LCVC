
function covidPT_OC

clear all; close all; clc; format long;

% ficheiro que resulta do NEOS - solucao OC 
M=xlsread('data-pt-covid19-15abril');

x=load('resultado_neos_umax025state060.txt');  % umax = 0.25 e state constraint <= 0.6

x2 = load('neos_umax05_state2_3.txt'); % umax = 0.5 e state constraint <= 2/3 Imax

t = x(:, 1); 
S = x(:, 2); 
A = x(:, 3); 
I = x(:, 4); 
P = x(:, 5); 
u = x(:, 6); % control

% ---- figures --------
N0 = 10295894 + 2 + 2/0.15;
time = 1:1:120; 

figure; 
hold on
plot(t, I,'m-', 'Linewidth', 2);
plot(x2(:, 1), x2(:, 4),'b--', 'LineWidth', 2); 
plot(time, M(time, 2)./N0,'r*', 'LineWidth', 2); 
title('Active infected I'); 
legend('u <= 0.25, 0.6*Imax', 'u <= 0.5, 2/3*Imax', 'real data'); 

figure; 
hold on; 
plot(u, 'm-', 'LineWidth', 2); 
plot(x2(:, 6), 'b--', 'LineWidth', 2); 
legend('u <= 0.25, 0.6*Imax', 'u <= 0.5, 2/3*Imax'); 
title('control'); 








