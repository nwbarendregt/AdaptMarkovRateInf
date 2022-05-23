% Figure_2_Generate.m
% Generates Figure 2 from Barendregt, et al., 2022.
clear
% Load performance data from adaptive and fixed-period algorithms for
% comparison:
load('tsb_adapt_sampled_rate_data.mat')
load('tsb_periodic_sampled_rate_data.mat')
% Plot average convergence times (Fig. 2B):
figure
plot(sample_period,ones(1,length(sample_period))*mean(N_samples_adapt_D),'k--')
hold on
plot(sample_period,mean(N_samples_periodic_D,2),'k-+','linewidth',5,'markersize',20)
xlim([min(sample_period),max(sample_period)])
% Plot average MSE for each transition rate (Fig. 2C):
figure
plot(sample_period,ones(1,length(sample_period))*mean(MSE_adapt_D(1,:)),'b--')
hold on
plot(sample_period,ones(1,length(sample_period))*mean(MSE_adapt_D(2,:)),'r--')
plot(sample_period,mean(MSE_periodic_D(:,:,1),2),'b-+','linewidth',5,'markersize',20)
plot(sample_period,mean(MSE_periodic_D(:,:,2),2),'r-+','linewidth',5,'markersize',20)
xlim([min(sample_period),max(sample_period)])
% Load performance data from adaptive algorithm for fixed true transition
% rates:
load('tsb_adapt_fixed_rate_data.mat')
% Plot average convergence time for fixed true transition rates (Fig. 2D):
figure
imagesc([h0_true(1) h0_true(end)],[h1_true(1) h1_true(end)],mean(N_samples_adapt_D,3)')
colormap gray; colorbar; set(gca,'ydir','normal')
% Plot average MSE for h0 inference (Fig. 2E):
figure
imagesc([h0_true(1) h0_true(end)],[h1_true(1) h1_true(end)],mean(MSE_adapt_D(:,:,:,1),3)')
colormap winter; colorbar; set(gca,'ydir','normal')
% Plot average MSE for h1 inference (Fig. 2F):
figure
imagesc([h0_true(1) h0_true(end)],[h1_true(1) h1_true(end)],mean(MSE_adapt_D(:,:,:,2),3)')
colormap autumn; colorbar; set(gca,'ydir','normal')