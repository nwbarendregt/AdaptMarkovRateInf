% bin_struct_fixed_size.m
% Calculates performance metrics for adaptive inference of arbitrary
% finite-state Markov chains with binary transition rates.
clear
% Set prior over transition rates and chain size:
p_link = 0.1:0.2:0.9; n_nodes = 3; max_edges = 6;
% Set covariance reduction tolerance and number of samples for each prior:
Cov_reduct = 1e-2; N_trials = 1e3;
% Construct transition rate mesh and indicies:
hi_mesh = [0 1]; h_ind = 1:max_edges;
% Construct infinitesimal transition matrix:
A = @(h) bin_struct_connect(h,n_nodes);
% Construct list of possible chain configurations:
states = dec2bin(0:2^(max_edges)-1)-'0';
% Pre-allocate performance matric storage:
N_samples = NaN(max_edges+1,length(p_link),N_trials);
L1_err = NaN(max_edges+1,length(p_link),N_trials);

for p = 1:length(p_link)
    % Construct prior over chain configurations:
    pn_0 = NaN(2^(max_edges),1);
    for i = 1:2^(max_edges)
        pn_0(i) = p_link(p)^sum(states(i,:))*(1-p_link(p))^(max_edges-sum(states(i,:)));
    end
    pn_0 = reshape(pn_0,2*ones(1,max_edges));
    for n_miss = 0:max_edges
        for n = 1:N_trials
            % Sample chain configuration given prior and number of missing
            % transitions:
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
            % Calculate convergence tolerance:
            V = det(Cov); Cov_tol = V*Cov_reduct;
            while V > Cov_tol
                Tp = T(end); Xp = X(end);
                % Calculate loose upper bound for optimization routine:
                ub = Tp+10/(sum(E(1+(Xp-1)*(n_nodes-1):1+(Xp-1)*(n_nodes-1)+(n_nodes-2))));
                % Calculate optimal time of next measurement given previous
                % measurements:
                Tn = fminbnd(@(T)var_bin_struct(T,A,pn,Xp,Tp,states,n_nodes,max_edges),Tp,ub);
                T = [T Tn];
                % Draw next measurement at optimal time:
                while Tn > sim_T
                    sim_Xp = sim_X;
                    [sim_X,sim_T] = bin_struct_sim(sim_X,sim_T,A(h_true));
                end
                X = [X sim_Xp];
                % Update posterior based on new measurement:
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
            N_samples(n_miss+1,p,n) = length(T);
            [~,I] = max(pn,[],'all','linear');
            h_pred = states(I,:);
            L1_err(n_miss+1,p,n) = sum(abs(h_pred-h_true))/max_edges;
        end
    end
end

save('bin_struct_fixed_size_data.m')