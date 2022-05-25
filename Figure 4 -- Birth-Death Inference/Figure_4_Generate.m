% Figure_4_Generate.m
% Generates Figure 4 from Barendregt et al., 2022.
clear
% Load performance metrics:
load('bd_adapt_fixed_rate_data.mat')
% Plot average convergence times for different pairs of transition rates
% (Fig. 4B):
figure
imagesc([mu_true(1) mu_true(end)],[lambda_true(1) lambda_true(end)],(mean(N_samples_adapt,3))',...
    'alphadata',~isnan(mean(N_samples_adapt,3,'omitnan')'))
colormap gray; colorbar; set(gca,'ydir','normal')
% Plot average MSE for inferring death rate mu (Fig. 4C):
figure
imagesc([mu_true(1) mu_true(end)],[lambda_true(1) lambda_true(end)],(mean(MSE_adapt(:,:,:,1),3))',...
    'alphadata',~isnan(mean(MSE_adapt(:,:,:,1),3,'omitnan')'))
colormap winter; colorbar; set(gca,'ydir','normal')
% Plot average MSE for inferring birth rate lambda (Fig. 4D):
figure
imagesc([mu_true(1) mu_true(end)],[lambda_true(1) lambda_true(end)],(mean(MSE_adapt(:,:,:,2),3))',...
    'alphadata',~isnan(mean(MSE_adapt(:,:,:,2),3,'omitnan')'))
colormap autumn; colorbar; set(gca,'ydir','normal')