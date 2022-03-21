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



