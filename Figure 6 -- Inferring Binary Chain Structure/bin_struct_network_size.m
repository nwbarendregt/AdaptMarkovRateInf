% bin_struct_network_size.m
% Calculates performance metrics for adaptive inference of arbitrary Markov
% chains with binary transition rates, fixed prior, and different state
% space sizes. Script is designed to test one state space size at a time.
clear
% Set prior over transition rates:
p_link = 0.5;
% Set state space size:
n_nodes = 3; max_edges = n_nodes*(n_nodes-1); ind = n_nodes-1;
% Set covariance reduction convergence condition and number of simulations
% for each state space size:
Cov_reduct = 1e-2; N_trials = 1e3;
% Construct transition rate mesh:
hi_mesh = [0 1]; h_ind = 1:max_edges;
% Construct infinitesimal generator matrix:
A = @(h) bin_struct_connect(h,n_nodes);
% Construct list of chain configurations:
states = dec2bin(0:2^(max_edges)-1)-'0';
% Pre-allocate performance metric storage:
N_samples = NaN(max_edges+1,N_trials);
L1_err = NaN(max_edges+1,N_trials);
% Construct prior over chain configurations based on transition rate prior:
pn_0 = NaN(2^(max_edges),1);
for i = 1:2^(max_edges)
    pn_0(i) = p_link^sum(states(i,:))*(1-p_link)^(max_edges-sum(states(i,:)));
end
pn_0 = reshape(pn_0,2*ones(1,max_edges));

for n = 1:N_trials
    N_samples_n = NaN(max_edges+1,1);
    L1_err_n = NaN(max_edges+1,1);
    for n_miss = 0:max_edges
        % Sample chain configuration based on state space size and prior:
        h_true = ones(1,max_edges); h_true(randperm(max_edges,n_miss))= 0;
        % Initialize simulation:
        X = randi(n_nodes); T = 0;
        sim_X = X; sim_T = T;
        % Calculate initial covariance determinant:
        pn = pn_0;
        E = NaN(1,max_edges); Cov = zeros(max_edges);
        for i = 1:max_edges
            if i ~= 2
                E(i) = sum(hi_mesh.*squeeze(sum(pn,h_ind(h_ind ~= i)))');
                Cov(i,i) = sum(hi_mesh.^2.*squeeze(sum(pn,h_ind(h_ind ~= i)))')-E(i)^2;
            else
                E(i) = sum(hi_mesh.*squeeze(sum(pn,h_ind(h_ind ~= i))));
                Cov(i,i) = sum(hi_mesh.^2.*squeeze(sum(pn,h_ind(h_ind ~= i))))-E(i)^2;
            end
        end
        % Calculate covariance convergence criteria:
        V = det(Cov); Cov_tol = V*Cov_reduct;
        while V > Cov_tol
            Tp = T(end); Xp = X(end);
            % Calculate loose upper bound for optimization routine:
            ub = Tp+10/(sum(E(1+(Xp-1)*(n_nodes-1):1+(Xp-1)*(n_nodes-1)+(n_nodes-2))));
            % Calculate optimal time to draw next measurement based on
            % previous measurements:
            Tn = fminbnd(@(T)var_bin_struct(T,A,pn,Xp,Tp,states,n_nodes,max_edges),Tp,ub);
            T = [T Tn];
            % Draw next measurement at optimal time:
            while Tn > sim_T
                sim_Xp = sim_X;
                [sim_X,sim_T] = bin_struct_sim(sim_X,sim_T,A(h_true));
            end
            X = [X sim_Xp];
            % Update posterior with new measurement:
            p_X = NaN(2^(max_edges),1);
            for i = 1:2^(max_edges)
                p_X(i) = bin_struct_trans_prob(A(states(i,:)),Tn-Tp,X(end),Xp);
            end
            p_X = reshape(p_X,2*ones(1,max_edges));
            pn = p_X.*pn; pn = pn/sum(pn,'all');
            % Calculate updated posterior covariance:
            E = NaN(1,max_edges); Cov = zeros(max_edges);
            for i = 1:max_edges
                if i ~= 2
                    E(i) = sum(hi_mesh.*squeeze(sum(pn,h_ind(h_ind ~= i)))');
                    Cov(i,i) = sum(hi_mesh.^2.*squeeze(sum(pn,h_ind(h_ind ~= i)))')-E(i)^2;
                else
                    E(i) = sum(hi_mesh.*squeeze(sum(pn,h_ind(h_ind ~= i))));
                    Cov(i,i) = sum(hi_mesh.^2.*squeeze(sum(pn,h_ind(h_ind ~= i))))-E(i)^2;
                end
            end
            V = det(Cov);
        end
        % Store performance metrics:
        N_samples_n(n_miss+1) = length(T);
        [~,I] = max(pn(:));
        h_pred = states(I,:);
        L1_err_n(n_miss+1) = sum(abs(h_pred-h_true))/max_edges;
    end
    N_samples(:,n) = N_samples_n;
    L1_err(:,n) = L1_err_n;
end
% Incorporate new data into data structure:
load('bin_struct_network_size.mat','bin_struct_data');
bin_struct_data(ind).n_nodes = n_nodes;
bin_struct_data_data(ind).N_samples = N_samples;
bin_struct_data(ind).L1_err = L1_err;
save('bin_struct_network_size.mat','bin_struct_data');