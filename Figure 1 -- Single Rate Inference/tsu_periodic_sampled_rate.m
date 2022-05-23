% tsu_periodic_sampled_rate.m
% Calculates performance metrics for fixed-period inference of single
% transition rates, averaged over sampled transition rates.
clear
% Load transition rate samples:
load('gamma_samples.mat');
% Define gamma prior used to generate transition rates:
prior = @(x) beta^alpha*x.^(alpha-1).*exp(-beta*x)/gamma(alpha);
% Set posterior variance convergence tolerance:
Var_tol = 0.1;
% Define integration mesh and prior:
h0_mesh = linspace(0,10,1000);
pn_0 = prior(h0_mesh);
% Define sample periods to test:
sample_period = linspace(0.1,1,10);
% Pre-allocate performance metric storage:
N_samples = NaN(length(sample_period),length(h0));
MSE = NaN(length(sample_period),length(h0));

for i = 1:length(h0)
    N_i = NaN(length(sample_period),1);
    MSE_i = NaN(length(sample_period),1);
    for j = 1:length(sample_period)
        Tn = sample_period(j);
        pn = pn_0;
        % Compute initial posterior covariance:
        E_h0 = trapz(h0_mesh,h0_mesh.*pn);
        V = trapz(h0_mesh,(h0_mesh-E_h0).^2.*pn);
        T = 0;
        while V >= Var_tol
            % Store time of next measurement given fixed period:
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
        N_i(j) = length(T);
        MSE_i(j) = trapz(h0_mesh,(h0_mesh-h0(i)).^2.*pn);
    end
    N_samples(:,i) = N_i;
    MSE(:,i) = MSE_i;
end

N_samples_periodic = N_samples;
MSE_periodic = MSE;

save('tsu_periodic_sampled_rate_data.mat');