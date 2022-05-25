% Figure_6_Generate.m
% Generates Figure 6 from Barendregt et al., 2022.
clear
% Load performance data for different transition rate priors:
load('bin_struct_fixed_size_data.mat')
% Plot average convergence speeds (Fig. 6B):
figure
plot(0:max_edges,flip(mean(N_samples,3)),'-*','linewidth',5,'markersize',20)
legend({'p=0.1','p=0.3','p=0.5','p=0.7','p=0.9'})
% Plot average L1 errors (Fig. 6C):
figure
plot(0:max_edges,flip(mean(L1_err,3)),'-*','linewidth',5,'markersize',20)
legend({'p=0.1','p=0.3','p=0.5','p=0.7','p=0.9'})
% Load performance data for different state space sizes:
load('bin_struct_network_size.mat')
f1 = figure; f2 = figure;
for i = 1:length(bin_struct_data)
    p_miss = linspace(0,1,i*(i+1)+1);
    N_samples = flip(mean(bin_struct_data(i).N_samples,2));
    % Plot average convergence speeds (Fig. 6D):
    figure(f1);
    plot(p_miss,N_samples,'-*','linewidth',5,'markersize',20); hold on
    err = flip(mean(bin_struct_data(i).L1_err,2));
    % Plot average L1 errors (Fig. 6E):
    figure(f2);
    plot(p_miss,err,'-*','linewidth',5,'markersize',20); hold on
end