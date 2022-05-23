% tsb_adapt_fixed_rate.m
% Calculates performance metrics for adaptive inference of two
% transition rates, averaged over realizations of each pair of transition
% rates.
clear
% Define bivariate gamma prior for transition rates:
m1 = 2; m2 = 2; a = 1; b = 1;
prior = @(x,y) bvgamma(x,y,a,b,m1,m2);
% Set posterior covariance determinant convergence tolerance and number of
% realizations:
Cov_tol = 0.1; N_trials = 1e3;
% Define true transition rate mesh:
h0_true = linspace(0,4,10); h1_true = linspace(0,4,10);
% Define integration mesh and prior:
h0_mesh = linspace(0,20); h1_mesh = linspace(0,20);
[H0_mesh,H1_mesh] = meshgrid(h0_mesh,h1_mesh);
pn_0 = prior(H0_mesh,H1_mesh); pn_0(1,1) = 0;
% Pre-allocate performance metric storage:
N_samples = NaN(length(h0_true),length(h1_true),N_trials);
MSE = NaN(length(h0_true),length(h1_true),N_trials,2);

for k = 1:N_trials
    N_k = NaN(length(h0_true),length(h1_true));
    MSE_k = NaN(length(h0_true),length(h1_true),2);
    for i = 1:length(h0_true)
        for j = 1:length(h1_true)
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
                    [sim_X,sim_T] = tsb_sim(sim_X,sim_T,h0_true(i),h1_true(j));
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
            N_k(i,j) = length(T);
            MSE_k(i,j,:) = [trapz(h0_mesh,((h0_mesh-h0_true(i)).^2).*trapz(h1_mesh,pn,1),2);...
                trapz(h1_mesh,((h1_mesh-h1_true(j)).^2)'.*trapz(h0_mesh,pn,2),1)];
        end
    end
    N_samples(:,:,k) = N_k;
    MSE(:,:,k,:) = MSE_k;
end

N_samples_adapt_D = N_samples;
MSE_adapt_D = MSE;

save('tsb_adapt_fixed_rate_data.mat');