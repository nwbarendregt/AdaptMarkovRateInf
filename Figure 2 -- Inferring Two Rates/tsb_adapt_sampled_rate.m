% tsb_adapt_sampled_rate.m
% Calculates performance metrics for adaptive inference of two
% transition rates, averaged over pairs of transition rates.
clear
% Load transition rate samples:
load('bvgamma_samples.mat');
% Define bivariate gamma prior for transition rates:
prior = @(x,y) bvgamma(x,y,a,b,m1,m2);
% Set posterior covariance determinant convergence tolerance:
Cov_tol = 0.1;
% Define integration mesh and prior:
h0_mesh = linspace(0,20); h1_mesh = linspace(0,20);
[H0_mesh,H1_mesh] = meshgrid(h0_mesh,h1_mesh);
pn_0 = prior(H0_mesh,H1_mesh); pn_0(1,1) = 0;
% Pre-allocate performance metric storage:
N_samples = NaN(1,length(h0));
MSE = NaN(2,length(h0));

for i = 1:length(h0)
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
        Tp = T(end); Xp = X(end);
        % Define loose upper bound to use for optimization:
        ub = Tp+10/(E_h0+E_h1);
        % Compute and store optimal time of next measurement:
        Tn = fminbnd(@(T)var_tsb(T,h0_mesh,h1_mesh,pn,Xp,Tp),Tp,ub);
        T = [T Tn];
        % Draw measurement at computed time:
        while Tn > sim_T
            [sim_X,sim_T] = tsb_sim(sim_X,sim_T,h0(i),h1(i));
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
    N_samples(i) = length(T);
    MSE(:,i) = [trapz(h0_mesh,((h0_mesh-h0(i)).^2).*trapz(h1_mesh,pn,1),2);...
        trapz(h1_mesh,((h1_mesh-h1(i)).^2)'.*trapz(h0_mesh,pn,2),1)];
end

N_samples_adapt_D = N_samples;
MSE_adapt_D = MSE;

save('tsb_adapt_sampled_rate_data.mat');