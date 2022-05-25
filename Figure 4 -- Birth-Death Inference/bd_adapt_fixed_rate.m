% bd_adapt_fixed_rate.m
% Calculates performance metrics for adaptive inference of a linear
% birth-death process, averaged over many realizations using the same fixed
% transition rates.
clear
% Set bivariate gamma prior for transition rates:
m1 = 2; m2 = 2; a = 1; b = 1;
prior = @(x,y) bvgamma(x,y,a,b,m1,m2);
% Set covariance convergence tolerance, number of trials for each pair of
% transition rates, and numerical truncation of infinite state space:
Cov_tol = 0.1; N_trials = 1e2; N_trunc = 10;
% Define true transition rates to test:
mu_true = linspace(eps,4,10); lambda_true = linspace(eps,4,10);
% Construct integration mesh and truncated prior:
mu_mesh = linspace(eps,20); lambda_mesh = linspace(eps,20);
[M_mesh,L_mesh] = meshgrid(mu_mesh,lambda_mesh);
pn_0 = 2*prior(M_mesh,L_mesh); pn_0(1,1) = 0;
for i = 1:length(mu_mesh)
    pn_0(i,1:i) = 0;
end
% Pre-allocate performance metric storage:
N_samples = NaN(length(mu_true),length(lambda_true),N_trials);
MSE = NaN(length(mu_true),length(lambda_true),N_trials,2);

for k = 1:N_trials
    N_k = NaN(1,length(lambda_true));
    MSE_k = NaN(1,length(lambda_true),2);
    for i = 1:length(mu_true)
        for j = 1:(i-1)
            % Initialize simulation:
            X = randi(N_trunc)-1; T = 0;
            sim_X = X; sim_T = T;
            % Calculate initial posterior covariance:
            pn = pn_0;
            E_mu = trapz(mu_mesh,trapz(lambda_mesh,M_mesh.*pn,1),2);
            E_lambda = trapz(lambda_mesh,trapz(mu_mesh,L_mesh.*pn,2),1);
            Cov = NaN(2,2);
            Cov(1,1) = trapz(mu_mesh,trapz(lambda_mesh,M_mesh.^2.*pn,1),2)-E_mu^2;
            Cov(1,2) = trapz(mu_mesh,trapz(lambda_mesh,M_mesh.*L_mesh.*pn,1),2)-E_mu*E_lambda; Cov(2,1) = Cov(1,2);
            Cov(2,2) = trapz(lambda_mesh,trapz(mu_mesh,L_mesh.^2.*pn,2),1)-E_lambda^2;
            V = det(Cov);
            while V >= Cov_tol
                Tp = T(end); Xp = X(end);
                % Compute loose upper bound for optimization:
                ub = Tp+10/(E_mu+E_lambda);
                % Calculate optimal time for next measurement given current
                % posterior:
                Tn = fminbnd(@(T)var_bd(T,mu_mesh,lambda_mesh,N_trunc,pn,Xp,Tp),Tp,ub);
                T = [T Tn];
                % Draw new measurement at optimal time:
                while Tn > sim_T
                    sim_Xp = sim_X;
                    [sim_X,sim_T] = bd_sim(sim_X,sim_T,mu_true(i),lambda_true(j));
                end
                X = [X sim_Xp];
                % Update posterior given new measurement:
                p = NaN(size(M_mesh));
                for ii = 1:length(mu_mesh)
                    for jj = 1:length(lambda_mesh)
                        p(ii,jj) = bd_trans_prob(mu_mesh(ii),lambda_mesh(jj),X(end),Xp,Tn,Tp);
                    end
                end
                p = p';
                pn = p.*pn; pn = pn/trapz(lambda_mesh,trapz(mu_mesh,pn,2),1);
                % Calculate updated posterior covariance:
                E_mu = trapz(mu_mesh,trapz(lambda_mesh,M_mesh.*pn,1),2);
                E_lambda = trapz(lambda_mesh,trapz(mu_mesh,L_mesh.*pn,2),1);
                Cov = NaN(2,2);
                Cov(1,1) = trapz(mu_mesh,trapz(lambda_mesh,M_mesh.^2.*pn,1),2)-E_mu^2;
                Cov(1,2) = trapz(mu_mesh,trapz(lambda_mesh,M_mesh.*L_mesh.*pn,1),2)-E_mu*E_lambda; Cov(2,1) = Cov(1,2);
                Cov(2,2) = trapz(lambda_mesh,trapz(mu_mesh,L_mesh.^2.*pn,2),1)-E_lambda^2;
                V = det(Cov);
            end
            % Store performance metrics:
            N_k(j) = length(T);
            MSE_k(:,j,:) = [trapz(mu_mesh,((mu_mesh-mu_true(i)).^2).*trapz(lambda_mesh,pn,1),2);...
                trapz(lambda_mesh,((lambda_mesh-lambda_true(j)).^2)'.*trapz(mu_mesh,pn,2),1)];
        end
    end
    N_samples(:,:,k) = N_k;
    MSE(:,:,k,:) = MSE_k;
end

N_samples_adapt = N_samples;
MSE_adapt = MSE;

save('bd_adapt_fixed_rate_data.mat','mu_true','lambda_true',...
    'N_samples_adapt','MSE_adapt')