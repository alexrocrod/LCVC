clear
close all
%
% USER INPUTS
%
L = 6; % linear system size
f_start = exp(1);
f_end = 1 + 1E-3;
max_rw_steps = 1E10; % maximum number of random walk steps per f value update
min_rw_steps = 1E4; % minimum random walk steps before checking flatness
p = 0.75; % flatness criteria
%
NN = 4; % number of nearest neighbours
N_atm = L^2; % total number of spins
%
M_list(:,1) = -N_atm : 2 : N_atm; % possible magnetization values
E_list(:,1) = 1/2 * N_atm * NN : -4 : -1/2 * N_atm * NN; % possible energy values
[ ~, index_M_max] = max(M_list); % to normalize the JDOS 
[ ~, index_E_min] = min(E_list); % to normalize the JDOS 
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
% preallocate variables
%
S_vector = nan(N_atm, 1);
H_EM = zeros(length(E_list), length(M_list));
JDOS = ones(length(E_list), length(M_list)); % JDOS values start as equal to one
log_JDOS = log(JDOS);
%
% main WL cycle
%
WL_timer = tic; % timer for WL sampling
%
% start in a random configuration
%
S_vector(:,1) = randsample([-1, 1], N_atm, true, [1/2, 1/2]);
%
f = f_start;
k = 1; % counter for rw steps
%
while k <= max_rw_steps && f >= f_end
    %
    M_old = sum(S_vector(:,1));
    E_old = function_Energy_Ising_2D_SS(L, S_vector);
    %
    % choose a random spin to flip
    %
    spin_flip_index = randperm(N_atm,1); % randomly choose spin to flip
    S_vector(spin_flip_index) = - S_vector(spin_flip_index); % flip spin
    %
    delta_E = - S_vector(spin_flip_index) .* ( ...
        S_vector(nnxpos(spin_flip_index)) + ...
        S_vector(nnxneg(spin_flip_index)) + ...
        S_vector(nnypos(spin_flip_index)) + ...
        S_vector(nnyneg(spin_flip_index))); % energy of bonds to NN
    %
    E_new = E_old + 2*delta_E; % calculate energy
    M_new = sum(S_vector); % calculate magnetization
    %
    if rand < min( [ exp(log_JDOS(E_list == E_old, M_list == M_old) - log_JDOS(E_list == E_new, M_list == M_new)), 1]) % accept
        %
        log_JDOS(E_list == E_new, M_list == M_new) = log(f) + log_JDOS(E_list == E_new, M_list == M_new); % update JDOS: multiply by f
        E_old = E_new;
        H_EM(E_list == E_new, M_list == M_new) = H_EM( E_list == E_new, M_list == M_new) + 1; % update histogram
        %
    else % reject
        %
        log_JDOS(E_list == E_old, M_list == M_old) = log(f) + log_JDOS(E_list == E_old, M_list == M_old); % update JDOS: multiply by f
        S_vector(spin_flip_index) = - S_vector(spin_flip_index); % unflip spin
        H_EM(E_list == E_old, M_list == M_old) = H_EM( E_list == E_old, M_list == M_old) + 1; % update histogram
        %
    end
    %
    converge_JDOS_max = max(H_EM(H_EM >0)) / mean(H_EM(H_EM >0)); % check flatness, using highest histogram value as reference
    converge_JDOS_min = min(H_EM(H_EM >0)) / mean(H_EM(H_EM >0)); % check flatness, using lowest histogram value as reference
    %
    if (k > min_rw_steps && converge_JDOS_max < (1 + (1-p)) && converge_JDOS_min > (1 - (1-p))) % histogram is flat
        %
        log_JDOS(log_JDOS == 0) = nan; % cleanup log JDOS of non-visited (E,M) values
        log_JDOS = log_JDOS - log_JDOS(index_E_min, index_M_max); % normalize JDOS in log scale
        JDOS = exp(log_JDOS); % JDOS in linear scale
        %
        H_EM = zeros(length(E_list), length(M_list)); % reset histogram
        disp(['histogram flatness reached for f = ',num2str(f),' in ',int2str(k),' steps'])
        %
        f = sqrt(f); % update f factor
        k = 0; % reset rw steps counter
        %
    end
    %
    k = k + 1;
    %
end
%
if k > max_rw_steps
    %
    disp(['convergence not reached for f = ',num2str(f),''])
    %
end
%
WL_time = toc(WL_timer); % register timer
disp(['WL time ', num2str(WL_time), ' seconds']); % display total WL time

%% plot3

hist_EM = JDOS;

figure(1)
hist_EM(hist_EM == 0) = nan;
b = bar3(log10(hist_EM));
zlabel('counts in log 10 scale')
% view(133,33) % angled view
view(2) % top-down view
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
% view(2)
cb = colorbar;
cb.Label.String="counts in log10 scale";


