% Figure_1_Generate.m
% Generates Figure 1 from Barendregt et al., 2022.
clear
% Load sampled performance metrics to compare adaptive inference to
% fixed-period inference.
load('tsu_adapt_sampled_rate_data.mat')
load('tsu_periodic_sampled_rate_data.mat')
% Plot average convergence time comparison (Fig. 1B):
figure
plot(sample_period,ones(1,length(sample_period))*mean(N_samples_adapt),'k--')
hold on
plot(sample_period,mean(N_samples_periodic,2),'k-+','linewidth',5,'markersize',20)
xlim([min(sample_period),max(sample_period)])
% Plot average MSE comparison (Fig. 1C):
figure
plot(sample_period,ones(1,length(sample_period))*mean(MSE_adapt),'k--')
hold on
plot(sample_period,mean(MSE_periodic,2),'k-+','linewidth',5,'markersize',20)
xlim([min(sample_period),max(sample_period)])
% Load adaptive inference performance metrics for fixed true transition
% rates:
load('tsu_adapt_fixed_rate_data.mat');
% Plot average convergence time (Fig. 1D):
figure
plot(h0_true,mean(N_samples_adapt,2),'k-*','linewidth',5,'markersize',20)
xlim([min(h0_true),max(h0_true)])
% Plot average MSE (Fig. 1E):
figure
plot(h0_true,mean(MSE_adapt,2),'k-*','linewidth',5,'markersize',20)
xlim([min(h0_true),max(h0_true)])