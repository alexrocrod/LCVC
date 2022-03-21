clear
close all
%
L = 4; % linear system size
REP = 1E5; % number of randomly generated spin configurations
%
NN = 4; % number of nearest neighbours
N_atm = L^2; % total number of spins
%
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = 1/2 * N_atm * NN : -4 : -1/2 * N_atm * NN; % possible energy values
%
% preallocate variables
%
S_vector = nan(N_atm, 1);
M = nan(REP, 1);
E = nan(REP, 1);
%
% main MC cycle
%
hist_E = nan(length(E_list), 2);
hist_E(:,1) = E_list;
%
hist_M = nan(length(M_list), 2);
hist_M(:,1) = M_list;
%
frames = [];
idx = 1;
for i= 1:9
    for j=1:log10(REP)
        frames(idx) = i*10^j;
        idx = idx+1;
    end
end
sort(frames)
%% 
% frames = [1 10 20 100 200 1000 2000 5000 1e4 1e5 1e6];
%
MC_timer = tic; % timer for MC sampling
%

fig = figure;
v = VideoWriter("video.avi");
open(v);

ax = gca;

for k = 1:REP
    %
    S_vector(:,1) = randi([-1, 0], N_atm, 1)';
    S_vector(S_vector(:,1) == 0, 1) = 1;
%     S_vector(:,1) = randsample([-1, 1], N_atm, true, [0.5,0.5]);
    M(k,1) = sum(S_vector(:,1));
    E(k,1) = function_Energy_Ising_2D_SS(L, S_vector);
    %
    idx = 1;
    if ismember(k,frames)
        disp(k)
%         figure(1)
        for E_index = 1:length(E_list)
            hist_E(E_index,2) = sum(E == E_list(E_index));
        end
        for M_index = 1:length(M_list)
            hist_M(M_index,2) = sum(M == M_list(M_index));
        end
        subplot(1,2,1)
        bar(hist_E(:,1),hist_E(:,2))
       
        subplot(1,2,2)
        bar(hist_M(:,1),hist_M(:,2))

        drawnow();
        F = getframe(ax);
        writeVideo(v,F);
%         clf(fig)

        idx = idx+1;
    end
end
% disp('movie:')
% movie(F);

%% save video


% for k = 1:length(F(:))
% %     movie(fig,F(k),1)
% %     frame = getframe(gcf);
%     writeVideo(v,F(k));
% end
close(gca)
close(v)


%% o

%
MC_time = toc(MC_timer); % register timer
disp(['MC time ', num2str(MC_time), ' seconds']); % display total MC time
%
% HISTOGRAMS


%
% (M,E) histogram
%
hist_timer = tic; % timer
%
hist_EM = nan(length(E_list), length(M_list));
%
for E_index = 1:length(E_list)
    %
    for M_index = 1:length(M_list)
        %
        hist_EM(E_index, M_index) = nnz(M == M_list(M_index) & E == E_list(E_index)) ;
        %
    end
    %
end
%
hist_time = toc(hist_timer); % register (E,M) histogram timer
disp(['(E,M) histogram time ', num2str(hist_time), ' seconds']); % display total histogram time
disp(['MC + histogram time ', num2str(MC_time + hist_time), ' seconds']); % display total time
%
% Histogram normalization to obtain the JDOS estimate
%
JDOS = hist_EM./REP * 2^N_atm;
%
% Histogram:
% 
% Magnetization
bar(hist_M(:,1),hist_M(:,2))
title('Magnetization histogram')
%
% Energy
figure(2)
bar(hist_E(:,1),hist_E(:,2))
title('Energy histogram')

getframe()


%% 3D
% joint E-M hist
figure(3)
contour(-E_list,M_list,hist_EM,100)
%
hist_EM(hist_EM==0)=nan;
colorbar
% view(2)
view(-45,45)
surf(-E_list,M_list,hist_EM)
%
colorbar 
view(-45,45)
% bar3(hist_EM)
% creat_bar3_plot()






