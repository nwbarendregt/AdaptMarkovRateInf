% tsu_adapt_sampled_rate.m
% Calculates performance metrics for adaptive inference of fixed true
% transition rates, averaged over realizations of the same rate.
clear
% Define gamma prior for transition rates:
alpha = 2; beta = 1;
prior = @(x) beta^alpha*x.^(alpha-1).*exp(-beta*x)/gamma(alpha);
% Set posterior variance convergence tolerance and number of realizations
% of each transition rate:
Var_tol = 0.1; N_trials = 1e3;
% Define mesh of true transition rates:
h0_true = linspace(0,4,10);
% Define integration mesh and prior:
h0_mesh = linspace(0,10,1000);
pn_0 = prior(h0_mesh);
% Pre-allocate performance metric storage:
N_samples = NaN(length(h0_true),N_trials);
MSE = NaN(length(h0_true),N_trials);

for j = 1:N_trials
    N_k = NaN(length(h0_true),1);
    MSE_k = NaN(length(h0_true),1);
    for i = 1:length(h0_true)
        pn = pn_0;
        % Compute initial posterior covariance:
        E_h0 = trapz(h0_mesh,h0_mesh.*pn);
        V = trapz(h0_mesh,(h0_mesh-E_h0).^2.*pn);
        T = 0;
        while V >= Var_tol
            % Define loose upper bound to use for optimization:
            ub = 10/E_h0;
            % Compute and store optimal time of next measurement:
            Tn = fminbnd(@(T)var_tsu(T,h0_mesh,pn),0,ub);
            T = [T T(end)+Tn];
            % Draw measurement at computed time and update posterior:
            if Tn < -log(rand)/h0_true(i)
                pn = exp(-h0_mesh*Tn).*pn;
                pn = pn/trapz(h0_mesh,pn);
            else
                pn = (1-exp(-h0_mesh*Tn)).*pn;
                pn = pn/trapz(h0_mesh,pn);
            end
            % Update posterior variance:
            E_h0 = trapz(h0_mesh,h0_mesh.*pn);
            V = trapz(h0_mesh,(h0_mesh-E_h0).^2.*pn);
        end
        % Calculate and store performance metrics:
        N_k(i) = length(T);
        MSE_k(i) = trapz(h0_mesh,(h0_mesh-h0_true(i)).^2.*pn);
    end
    N_samples(:,j) = N_k;
    MSE(:,j) = MSE_k;
end

N_samples_adapt = N_samples;
MSE_adapt = MSE;

save('tsu_adapt_fixed_rate_data.mat');