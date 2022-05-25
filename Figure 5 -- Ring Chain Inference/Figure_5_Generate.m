% Figure_5_Generate.m
% Generates Figure 5 from Barendregt et al., 2022.
clear
% Load performance metrics for different state space sizes:
load('ring_sampled_rate.mat')
% Pre-allocate performance metric storage:
N_samples = NaN(1,length(sampled_ring_data)); MSE = NaN(2,length(sampled_ring_data));
for i = 1:length(sampled_ring_data)
    N_samples(i)=mean(sampled_ring_data(i).N_samples);
    MSE(:,i) = mean(sampled_ring_data(i).MSE,2);
end
% Plot average convergence times for different state space sizes (Fig. 5B):
figure
plot((1:length(sampled_ring_data))+2,N_samples,'k-*','linewidth',5,'markersize',20)
% Plot average MSE for different state space sizes (Fig. 5C):
figure
plot((1:length(sampled_ring_data))+2,MSE(1,:),'b-*','linewidth',5,'markersize',20)
hold on
plot((1:length(sampled_ring_data))+2,MSE(2,:),'r-*','linewidth',5,'markersize',20)
% Load performance metrics for fixed state space size, averaged over fixed
% transition rates:
load('ring_adapt_fixed_rate_data.mat');
% Plot average convergence times for different fixed pairs of transition
% rates (Fig. 5D):
figure
imagesc([hp_true(1) hp_true(end)],[hm_true(1) hm_true(end)],mean(N_samples_adapt,3)')
colormap gray; colorbar; set(gca,'ydir','normal')
% Plot average MSE associated with inferring h_+ (Fig. 5E):
figure
imagesc([hp_true(1) hp_true(end)],[hm_true(1) hm_true(end)],mean(MSE_adapt(:,:,:,1),3)')
colormap winter; colorbar; set(gca,'ydir','normal')
% Plot average MSE associated with inferring h_-:
figure
imagesc([hp_true(1) hp_true(end)],[hm_true(1) hm_true(end)],mean(MSE_adapt(:,:,:,2),3)')
colormap autumn; colorbar; set(gca,'ydir','normal')