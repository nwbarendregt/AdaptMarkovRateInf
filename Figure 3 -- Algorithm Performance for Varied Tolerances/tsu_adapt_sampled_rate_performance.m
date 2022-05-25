% tsu_adapt_sampled_rate.m
% Calculates performance metrics for adaptive inference of single
% transition rates, averaged over sampled transition rates, for different
% convergence tolerances.
clear
% Load transition rate samples:
load('gamma_samples.mat');
% Define gamma prior used to generate transition rates:
prior = @(x) beta^alpha*x.^(alpha-1).*exp(-beta*x)/gamma(alpha);
% Set posterior variance convergence tolerance:
Var_tol = [0.1 0.05 0.01 0.005 0.001];
% Define integration mesh and prior:
h0_mesh = linspace(0,10);
pn_0 = prior(h0_mesh);
% Pre-allocate performance metric storage:
N_samples = NaN(length(h0),length(Var_tol));
MSE = NaN(length(h0),length(Var_tol));

for i = 1:length(h0)
    N_samples_i = NaN(1,length(Var_tol));
    MSE_i = NaN(1,length(Var_tol));
    for j = 1:length(Var_tol)
        % Compute initial posterior variance:
        pn = pn_0;
        E_h0 = trapz(h0_mesh,h0_mesh.*pn);
        V = trapz(h0_mesh,(h0_mesh-E_h0).^2.*pn);
        T = 0;
        while V >= Var_tol(j)
            % Define loose upper bound to use for optimization:
            ub = 10/E_h0;
            % Compute and store optimal time of next measurement:
            Tn = fminbnd(@(T)var_tsu(T,h0_mesh,pn),0,ub);
            T = [T T(end)+Tn];
            % Draw measurement at computed time and update posterior:
            if Tn < -log(rand)/h0(i)
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
        N_samples_i(j) = length(T);
        MSE_i(j) = trapz(h0_mesh,(h0_mesh-h0(i)).^2.*pn);
    end
    N_samples(i,:) = N_samples_i;
    MSE(i,:) = MSE_i;
end

N_samples_adapt = N_samples;
MSE_adapt = MSE;

save('tsu_adapt_sampled_rate_performance_data.mat');