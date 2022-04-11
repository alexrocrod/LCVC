

%% sec1
close
clear


%% sec2
load res.txt
plot(res(:,1),res(:,2))
xlabel("Iterações Monte Carlo")
ylabel("Energia")