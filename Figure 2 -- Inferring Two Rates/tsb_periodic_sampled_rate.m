% tsb_periodic_sampled_rate.m
% Calculates performance metrics for fixed-period inference of two
% transition rates, averaged over pairs of transition rates.
clear
% Load transition rate samples:
load('bvgamma_samples.mat');
% Define bivariate gamma prior for transition rates:
prior = @(x,y) bvgamma(x,y,a,b,m1,m2);
% Set posterior covariance determinant convergence tolerance:
Cov_tol = 0.1;
% Define integration mesh and prior:
h0_mesh = linspace(0,20,1000); h1_mesh = linspace(0,20,1000);
[H0_mesh,H1_mesh] = meshgrid(h0_mesh,h1_mesh);
pn_0 = prior(H0_mesh,H1_mesh); pn_0(1,1) = 0;
% Define sampling periods to test:
sample_period = linspace(0.1,2,10);
% Pre-allocate performance metric storage:
N_samples = NaN(length(sample_period),length(h0));
MSE = NaN(length(sample_period),length(h0),2);

for i = 1:length(h0)
    N_i = NaN(length(sample_period),1);
    MSE_i = NaN(length(sample_period),2);
    for j = 1:length(sample_period)
        % Set sampling period:
        dT = sample_period(j);
        % Initialize simulation:
        X = 0; T = 0; sim_T = 0; sim_X = 0;
        % Compute initial posterior covariance and determinant:
        pn = pn_0;
        E_h0 = trapz(h0_mesh,trapz(h1_mesh,H0_mesh.*pn,1),2);
        E_h1 = trapz(h1_mesh,trapz(h0_mesh,H1_mesh.*pn,2),1);
        Cov = NaN(2,2);
        Cov(1,1) = trapz(h0_mesh,trapz(h1_mesh,H0_mesh.^2.*pn,1),2)-E_h0^2;
        Cov(1,2) = trapz(h0_mesh,trapz(h1_mesh,H0_mesh.*H1_mesh.*pn,1),2)-E_h0*E_h1; Cov(2,1) = Cov(1,2);
        Cov(2,2) = trapz(h1_mesh,trapz(h0_mesh,H1_mesh.^2.*pn,2),1)-E_h1^2;
        V = det(Cov);
        while V >= Cov_tol
            % Compute and store time of next measurement given sampling 
            % period:
            Tn = T(end)+dT;
            T = [T Tn];
            % Draw measurement at computed time:
            while Tn > sim_T
                [sim_X,sim_T] = tsb_sim(sim_X,sim_T,h0(j),h1(j));
            end
            X = [X 1-sim_X];
            % Update posterior given new measurement:
            if X(end) == 0
                if X(end-1) == 0
                    p = H1_mesh./(H0_mesh+H1_mesh)+H0_mesh.*exp(-(H0_mesh+H1_mesh)*(T(end)-T(end-1)))./(H0_mesh+H1_mesh); p(1,1) = 1;
                    pn = p.*pn; pn = pn/trapz(h1_mesh,trapz(h0_mesh,pn,2),1);
                else
                    p = H0_mesh./(H0_mesh+H1_mesh)-H0_mesh.*exp(-(H0_mesh+H1_mesh)*(T(end)-T(end-1)))./(H0_mesh+H1_mesh); p(1,1) = 0;
                    pn = p.*pn; pn = pn/trapz(h1_mesh,trapz(h0_mesh,pn,2),1);
                end
            else
                if X(end-1) == 0
                    p = H1_mesh./(H0_mesh+H1_mesh)-H1_mesh.*exp(-(H0_mesh+H1_mesh)*(T(end)-T(end-1)))./(H0_mesh+H1_mesh); p(1,1) = 0;
                    pn = p.*pn; pn = pn/trapz(h1_mesh,trapz(h0_mesh,pn,2),1);
                else
                    p = H0_mesh./(H0_mesh+H1_mesh)+H1_mesh.*exp(-(H0_mesh+H1_mesh)*(T(end)-T(end-1)))./(H0_mesh+H1_mesh); p(1,1) = 1;
                    pn = p.*pn; pn = pn/trapz(h1_mesh,trapz(h0_mesh,pn,2),1);
                end
            end
            % Update posterior variance and determinant:
            E_h0 = trapz(h0_mesh,trapz(h1_mesh,H0_mesh.*pn,1),2);
            E_h1 = trapz(h1_mesh,trapz(h0_mesh,H1_mesh.*pn,2),1);
            Cov = NaN(2,2);
            Cov(1,1) = trapz(h0_mesh,trapz(h1_mesh,H0_mesh.^2.*pn,1),2)-E_h0^2;
            Cov(1,2) = trapz(h0_mesh,trapz(h1_mesh,H0_mesh.*H1_mesh.*pn,1),2)-E_h0*E_h1; Cov(2,1) = Cov(1,2);
            Cov(2,2) = trapz(h1_mesh,trapz(h0_mesh,H1_mesh.^2.*pn,2),1)-E_h1^2;
            V = det(Cov);
        end
        % Calculate and store performance metrics:
        N_i(j) = length(T);
        MSE_i(j,:) = [trapz(h0_mesh,((h0_mesh-h0(i)).^2).*trapz(h1_mesh,pn,1),2);...
            trapz(h1_mesh,((h1_mesh-h1(i)).^2)'.*trapz(h0_mesh,pn,2),1)];
    end
    N_samples(:,i) = N_i;
    MSE(:,i,:) = MSE_i;
end

N_samples_periodic_D = N_samples;
MSE_periodic_D = MSE;

save('tsb_periodic_sampled_rate_data.mat');