% ring_adapt_sampled_rate_network_size.m
% Calculates performance metrics for adaptive inference of a ring chain of
% states for different state space sizes. Script is designed to test one
% size at a time.
clear
% Load sampled transition rates for testing:
load('bvgamma_samples.mat');
hp = hp(1:100); hm = hm(1:100);
% Set state space size to test:
m = 3; ind = m-2;
% Construct infinitesimal generator matrix:
A = @(x,y) diag(-(x+y)*ones(1,m))+diag(x*ones(1,m-1),1)+...
    diag(y*ones(1,m-1),-1)+x*([zeros(m-1,1);1]*[1 zeros(1,m-1)])+...
    y*([1;zeros(m-1,1)]*[zeros(1,m-1) 1]);
% Set bivariate gamma prior for transition rates:
prior = @(x,y) bvgamma(x,y,a,b,m1,m2);
% Set convergence tolerance:
Cov_tol = 0.1;
% Construct integration mesh and prior:
hp_mesh = linspace(0,20); hm_mesh = linspace(0,20);
[Hp_mesh,Hm_mesh] = meshgrid(hp_mesh,hm_mesh);
pn_0 = prior(Hp_mesh,Hm_mesh); pn_0(1,1) = 0;
% Pre-allocate performance metric storage:
N_samples = NaN(1,length(hp));
MSE = NaN(2,length(hp));

for i = 1:length(hp)
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
        % Calculate loose upper bound for optimization routine:
        ub = Tp+10/(E_hp+E_hm);
        % Determine optimal time to draw next measurement based on previous
        % measurements:
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
    N_samples(i) = length(T);
    MSE(:,i) = [trapz(hp_mesh,((hp_mesh-hp(i)).^2).*trapz(hm_mesh,pn,1),2) ...
        trapz(hm_mesh,((hm_mesh-hm(i)).^2)'.*trapz(hp_mesh,pn,2),1)];
end
% Incorporate new data into data structure:
% load('ring_sampled_rate.mat','sampled_ring_data');
% sampled_ring_data(ind).N_states = N_states;
% sampled_ring_data(ind).N_samples = N_samples;
% sampled_ring_data(ind).MSE = MSE;
save('ring_sampled_rate.mat','sampled_ring_data');