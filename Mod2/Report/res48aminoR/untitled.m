close all
clear all
clc


load res.txt

plot(res(:,1),res(:,2))
xlabel("Iterações Monte Carlo")
ylabel("Energia")