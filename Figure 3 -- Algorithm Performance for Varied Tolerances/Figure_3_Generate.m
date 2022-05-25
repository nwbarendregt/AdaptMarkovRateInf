% Figure_3_Generate.m
% Generates Figure 3 from Barendregt et al., 2022.
clear
% Load performance metrics for inference on simple networks:
load('tsu_adapt_sampled_rate_performance_data.mat')
load('tsb_adapt_sampled_rate_performance_data.mat')
% Plot average convergence time for different convergence tolerances (Fig.
% 3A):
figure
loglog(Var_tol,mean(N_samples_adapt,1),'k-*','linewidth',5,'markersize',20) 
hold on
loglog(Cov_tol,mean(N_samples_adapt_D,1),'k-s','linewidth',5,'markersize',20,...
    'markerfacecolor','k')
% Plot average MSE for different convergence tolerances (Fig. 3B):
figure
loglog(Var_tol,mean(MSE_adapt,1),'k-*','linewidth',5,'markersize',20)
hold on
MSE_adapt_D = squeeze(mean(MSE_adapt_D,1));
loglog(Cov_tol,MSE_adapt_D(1,:),'b-s','linewidth',5,'markersize',20,...
    'markerfacecolor','b')
loglog(Cov_tol,MSE_adapt_D(2,:),'r-s','linewidth',5,'markersize',20,...
    'markerfacecolor','r')