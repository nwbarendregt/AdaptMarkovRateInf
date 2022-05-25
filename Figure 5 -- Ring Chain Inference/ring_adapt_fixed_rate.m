% ring_adapt_fixed_rate.m
% Calculates performance metrics for adaptive inference on a ring chain
% of states with clockwise (h_+) and counter-clockwise (h_-) transition
% rates.
clear
% Set bivariate gamma prior for transition rates:
m1 = 2; m2 = 2; a = 1; b = 1;
prior = @(x,y) bvgamma(x,y,a,b,m1,m2);
% Set convergence tolerance and number of samples for each pair of
% transition rates:
Cov_tol = 0.1; N_trials = 1e2;
% Set true transition rate mesh to test:
hp_true = linspace(0,4,10); hm_true = linspace(0,4,10);
% Define integration mesh and prior:
hp_mesh = linspace(0,20); hm_mesh = linspace(0,20);
[Hp_mesh,Hm_mesh] = meshgrid(hp_mesh,hm_mesh);
pn_0 = prior(Hp_mesh,Hm_mesh); pn_0(1,1) = 0;
% Set state space size and construct infinitesimal generator matrix:
m = 8;
A = @(x,y) diag(-(x+y)*ones(1,m))+diag(x*ones(1,m-1),1)+...
    diag(y*ones(1,m-1),-1)+x*([zeros(m-1,1);1]*[1 zeros(1,m-1)])+...
    y*([1;zeros(m-1,1)]*[zeros(1,m-1) 1]);
% Pre-allocate performance metric storage:
N_samples = NaN(length(hp_true),length(hm_true),N_trials);
MSE = NaN(length(hp_true),length(hm_true),N_trials,2);

for k = 1:N_trials
    N_k = NaN(length(hp_true),length(hm_true));
    MSE_k = NaN(length(hp_true),length(hm_true),2);
    for i = 1:length(hp_true)
        for j = 1:length(hm_true)
            % Initialize simulation:
            X = randi(m); T = 0;
            sim_X = X; sim_T = T;
            % Calculate initial covariance determinant:
            pn = pn_0;
            E_hp = trapz(hp_mesh,trapz(hm_mesh,Hp_mesh.*pn,1),2);
            E_hm = trapz(hm_mesh,trapz(hp_mesh,Hm_mesh.*pn,2),1);
            Cov = NaN(2,2);
            Cov(1,1) = trapz(hp_mesh,trapz(hm_mesh,Hp_mesh.^2.*pn,1),2)-E_hp^2;
            Cov(1,2) = trapz(hp_mesh,trapz(hm_mesh,Hp_mesh.*Hm_mesh.*pn,1),2)-E_hp*E_hm; Cov(2,1) = Cov(1,2);
            Cov(2,2) = trapz(hm_mesh,trapz(hp_mesh,Hm_mesh.^2.*pn,2),1)-E_hm^2;
            V = det(Cov);
            while V >= Cov_tol
                Tp = T(end); Xp = X(end);
                % Set loose upper bound for optimization routine:
                ub = Tp+10/(E_hp+E_hm);
                % Calculate optimal time to draw next measurement based on
                % previous measuerments:
                Tn = fminbnd(@(T)var_ring(T,hp_mesh,hm_mesh,A,m,pn,Xp,Tp),Tp,ub);
                T = [T Tn];
                % Draw next measurement at optimal time:
                while Tn > sim_T
                    sim_Xp = sim_X;
                    [sim_X,sim_T] = ring_sim(sim_X,sim_T,hp(i),hm(i),m);
                end
                X = [X sim_Xp];
                % Update posterior based on new measurement:
                p = NaN(size(Hp_mesh));
                for ii = 1:length(hp_mesh)
                    for jj = 1:length(hm_mesh)
                        p_h = ring_trans_prob(A(hp_mesh(ii),hm_mesh(jj)),Tn-Tp,Xp);
                        p(ii,jj) = p_h(X(end));
                    end
                end
                p = p';
                pn = p.*pn; pn = pn/trapz(hm_mesh,trapz(hp_mesh,pn,2),1);
                % Calculate updated posterior covariance:
                E_hp = trapz(hp_mesh,trapz(hm_mesh,Hp_mesh.*pn,1),2);
                E_hm = trapz(hm_mesh,trapz(hp_mesh,Hm_mesh.*pn,2),1);
                Cov = NaN(2,2);
                Cov(1,1) = trapz(hp_mesh,trapz(hm_mesh,Hp_mesh.^2.*pn,1),2)-E_hp^2;
                Cov(1,2) = trapz(hp_mesh,trapz(hm_mesh,Hp_mesh.*Hm_mesh.*pn,1),2)-E_hp*E_hm; Cov(2,1) = Cov(1,2);
                Cov(2,2) = trapz(hm_mesh,trapz(hp_mesh,Hm_mesh.^2.*pn,2),1)-E_hm^2;
                V = det(Cov);
            end
            % Store performance metrics:
            N_k(i,j) = length(T);
            MSE_k(i,j,:) = [trapz(hp_mesh,((hp_mesh-hp_true(i)).^2).*trapz(hm_mesh,pn,1),2);...
                trapz(hm_mesh,((hm_mesh-hm_true(j)).^2)'.*trapz(hp_mesh,pn,2),1)];
        end
    end
    N_samples(:,:,k) = N_k;
    MSE(:,:,k,:) = MSE_k;
end

N_samples_adapt = N_samples;
MSE_adapt = MSE;

save('ring_adapt_fixed_rate_data.mat');