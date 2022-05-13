clc
clear
close all

%% Read results from file

fileID=fopen('UcontrolVPDF2.txt','r');
formatSpec='%f %f %f %f %f %f %f';
xsize=[7 Inf];
x=fscanf(fileID,formatSpec,xsize);
x=x';

t = x(:, 1); 
s = x(:, 2); 
p = x(:, 3); 
a = x(:, 4); 
h = x(:, 5);
r = x(:, 6); 
u = x(:, 7);  % control

%% Making sure S+P+A+H+R=1
erro1 = mean(sum(x(:,2:6),2) - ones(size(x(:,2:6),1),1));
fprintf("Erro:%d\n",erro1)

%% Plots
figure; 
plot(t, s, t, p, t, a, t, h, t, r, 'LineWidth',1.5); 
legend('s', 'p', 'a', 'h', 'r'); 

figure; 
plot(t, p, t, a, t, h, t, r, 'LineWidth',1.5); 
legend('p', 'a', 'h', 'r');

figure; 
plot(t, a, t, h, t, r, 'LineWidth',1.5); 
legend('a', 'h', 'r');

figure; 
plot(t, u, 'LineWidth',1.5); 
legend('u');

