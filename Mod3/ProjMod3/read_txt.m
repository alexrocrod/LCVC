clc
clear
close all

%% Read results from file

addpath("v0\","VariaParam\Matlab\","VariaParam\Ucontrol\" ,"VariaInicial\Matlab\","VariaInicial\Ucontrol\")

% name = "PDF2"; % parametros e populacao do artigo %ParamNorm é igual
% name = "PDF3"; % parametros do artigo com situação inicial má
% name = "PDF4"; % parametros do artigo com situação inicial boa
% name = "ParamBom"; % situacao inicial do artigo com bons parametros
name = "ParamMau"; % situacao inicial do artigo com maus parametros


if (name=="PDF2" || name=="PDF3" || name=="PDF4")
    fig1 = sprintf("VariaInicial/Ucontrol/Geral%s.jpg",name);
    fig2 = sprintf("VariaInicial/Ucontrol/I2%s.jpg",name);
    fig3 = sprintf("VariaInicial/Ucontrol/%s.jpg",name);
elseif (name=="ParamBom" || name=="ParamMau" || name=="ParamNorm")
    fig1 = sprintf("VariaParam/Ucontrol/Geral%s.jpg",name);
    fig2 = sprintf("VariaParam/Ucontrol/I2%s.jpg",name);
    fig3 = sprintf("VariaParam/Ucontrol/%s.jpg",name);
else
    fprintf("Nome inválido")
    return
end


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
figure;
plot(t, s, t, p, t, a, t, h, t, r, 'LineWidth',1.5); 
legend('s', 'p', 'a', 'h', 'r'); 
saveas(gcf,fig1)


figure;
plot(t, p, t, a, t, h, t, r, 'LineWidth',1.5); 
legend('p', 'a', 'h', 'r');
saveas(gcf,fig2)

%% Plots
figure; 

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



