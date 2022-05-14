clc
clear
close all

%% Read results from file
% name = "PDF4";
% fig1 = sprintf("VariaInicial/Ucontrol/Geral%s.jpg",name);
% fig2 = sprintf("VariaInicial/Ucontrol/I2%s.jpg",name);
% fig3 = sprintf("VariaInicial/Ucontrol/%s.jpg",name);

name = "ParamBom";
fig1 = sprintf("VariaParam/Ucontrol/Geral%s.jpg",name);
fig2 = sprintf("VariaParam/Ucontrol/I2%s.jpg",name);
fig3 = sprintf("VariaParam/Ucontrol/%s.jpg",name);


fileID=fopen(sprintf('UcontrolV%s.txt',name),'r');
formatSpec='%f %f %f %f %f %f %f';
xsize=[7 Inf];
x=fscanf(fileID,formatSpec,xsize);
x=x';
fclose(fileID);

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

%% Plots init
fh = figure();
% fh.WindowState = 'maximized';
plot(t, s, t, p, t, a, t, h, t, r, 'LineWidth',1.5); 
legend('s', 'p', 'a', 'h', 'r'); 
saveas(gcf,fig1)
% 
fh = figure();
% fh.WindowState = 'maximized'; 
plot(t, p, t, a, t, h, t, r, 'LineWidth',1.5); 
legend('p', 'a', 'h', 'r');
saveas(gcf,fig2)
% 
% figure; 
% plot(t, a, t, h, t, r, 'LineWidth',1.5); 
% legend('a', 'h', 'r');
% 
% figure; 
% plot(t, u, 'LineWidth',1.5); 
% legend('u');


%% Plots
fh = figure();
fh.WindowState = 'maximized'; 
cols = 3;
rows = 2;

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

subplot(rows,cols,6)
plot(t,u,'LineWidth',1.5); 
title("Controlo (u)") 

saveas(gcf,fig3)



