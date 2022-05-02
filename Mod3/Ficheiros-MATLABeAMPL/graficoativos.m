function graficoativos

clear all; close all; clc; 

%  -- read excel fite with data ----
M=xlsread('data-pt-covid19-15abril');

tf = length(M);
time=1:1:tf;

ativos = M(1:tf,2)

% --- figure -----
figure; 
hold on
plot(time, ativos, 'r*'); 
x=xlabel('Time (days)');
y=ylabel('Active infected individuals - COVID-19');
t=title('COVID-19 Portuguese data: 02-03-202 -- 14-05-2021', 'FontSize', 12); 
set(y, 'FontName', 'Arial', 'FontSize', 12);
set(x, 'FontName', 'Arial', 'FontSize', 12);
set(t, 'FontName', 'Arial', 'FontSize', 12);
