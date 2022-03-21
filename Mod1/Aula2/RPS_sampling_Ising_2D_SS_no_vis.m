clear
close all
%
L = 8; % linear size of system
REP = 1E6; % number of RPS sweeps // equivalent to (N_atm + 1)*REP spin configurations
%
NN = 4; % number of nearest neighbours
N_atm = L^2; % total number of spins
%
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = 1/2 * N_atm * NN : -4 : -1/2 * N_atm * NN; % possible energy values
%
% PRELIMINARY CALCULATION OF NEIGHBOUR TABLES
%
nnxpos = nan(L^2,1); % right neighbout in x
nnxneg = nan(L^2,1); % left neighbour in x
nnypos = nan(L^2,1); % up neighbour in y
nnyneg = nan(L^2,1); % down neighbour in z 
%
for i=1:L % loop through all row positions
    %
    for j=1:L % loop through all column positions
        %
        [nnxpos(j+(i-1)*L), nnxneg(j+(i-1)*L), nnypos(j+(i-1)*L), nnyneg(j+(i-1)*L)] = function_NN_list_2D_SS(L,i,j);
        %
    end
    %
end
%
E_all = nan(REP, length(M_list)); % full Energy matrix of all sweeps
E_all(:,1) = -1/2 * N_atm * NN; % energy of all spins pointing up
E_all(:,length(M_list)) = -1/2 * N_atm * NN; % energy of all spins pointing down
%
% RPS SWEEPS
%
RPS_timer = tic; % timer for RPS sampling
%
for k = 1:REP % loop through all requested RPS loops
    %
    S_vector = ones(N_atm, 1); % vector with spins
    SFV(:,1) = randperm(N_atm); % spin flip vector (sequence of spins to flip)
    %
    for q = 2:N_atm % loop through magnetization values
        %
        S_vector(SFV(q-1)) = -1; % flip the spin
        %
        E_new = - S_vector(SFV(q-1)) .* ( ...
            S_vector(nnxpos((SFV(q-1)))) + ...
            S_vector(nnxneg((SFV(q-1)))) + ...
            S_vector(nnypos((SFV(q-1)))) + ...
            S_vector(nnyneg((SFV(q-1))))); % energy of bonds to NN
        %
        E_all(k, q) = E_all(k, q-1) + 2*E_new; % build the energy matrix
        %
    end
    %
end
%
RPS_time = toc(RPS_timer); % register timer
disp(['RPS time ', num2str(RPS_time), ' seconds']); % display total RPS time
%
% energy histogram
%
hist_E = nan(numel(E_list), 2);
hist_E(:,1) = E_list;
%
for E_index = 1:numel(E_list)
    %
    hist_E(E_index,2) = sum(sum(E_all == E_list(E_index)));
    %
end
%
% magnetization histogram 
%
hist_M = nan(numel(M_list), 2);
hist_M(:,1) = M_list;
hist_M(:,2) = REP;
%
% (M,E) histogram
%
hist_timer = tic; % timer
%
hist_EM = nan(numel(E_list), numel(M_list));
%
for E_index = 1:numel(E_list)
    %
    for M_index = 1:numel(M_list)
        %
        hist_EM(E_index, M_index) = nnz(E_all(:, M_index) == E_list(E_index)) ;
        %
    end
    %
end
%
hist_time = toc(hist_timer); % register (E,M) histogram timer
disp(['(E,M) histogram time ', num2str(hist_time), ' seconds']); % display total histogram time
disp(['RPS + histogram time ', num2str(RPS_time + hist_time), ' seconds']); % display total time
%
% Histogram normalization to obtain the JDOS estimate
%
eval(['load norm_factor_Ising_Natm_',int2str(N_atm),'.mat'])
%
JDOS = nan(numel(E_list), numel(M_list)); 
%
for q = 1:(N_atm+1)
    %
    JDOS(:,q) = hist_EM(:,q)/REP * norm_factor(q);
    %
end


figure(1)
hist_EM(hist_EM == 0) = nan;
b = bar3(log10(hist_EM));
zlabel('counts in log 10 scale')
view(133,33) % angled view
% view(2) % top-down view
axis([0 length(M_list)+1 0 length(E_list)+1])
xticks([1 (N_atm/4+1) (N_atm/2+1) (3*N_atm/4+1) N_atm+1])
yticks([1 ((length(E_list)-1)/4 + 1) ((length(E_list)-1)/2+1) (3*(length(E_list)-1)/4+1) length(E_list)])
%
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
%
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'XTickLabel', linspace(M_list(min(xt)), M_list(max(xt)), numel(xt))), xlabel('M')
set(gca, 'YTickLabel', linspace(E_list(min(yt)), E_list(max(yt)), numel(yt))), ylabel('E')
%
h = get(gca,'children');
%
for w = 1:length(h)
    %
    hc = get(h(w),'cdata');
    hz = get(h(w),'zdata');
    %
    for u = 1:(length(hc(:,1))/6)
        %
        if sum(nansum( hc(1+(u-1)*6 : u*6, 1:4) )) == 0
            %
            hc_new = hc;
            hc_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'cdata', hc_new);
            hc = hc_new;
            %
            hz_new = hz;
            hz_new(1+(u-1)*6 : u*6, 1:4) = nan(6,4);
            set(h(w), 'zdata', hz_new);
            hz = hz_new;
            %
        end
        %
    end
    %
end
% 
view(2)
cb = colorbar;
cb.Label.String="counts in log10 scale";
