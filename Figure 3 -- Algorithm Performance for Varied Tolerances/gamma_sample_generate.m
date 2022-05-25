% gamma_sample_generate.m
% Generates sample of transition rates from gamma prior to compare adaptive
% inference performance to fixed-period inference.
clear
% Define distribution parameters:
alpha = 2; beta = 1;
% Set number of transition rate samples:
N_Samples = 1e3;
% Generate samples:
h0 = sum(-log(rand(alpha,N_Samples)),1)/beta;
% Save transition rates for performance comparisons:
save('gamma_samples.mat','alpha','beta','h0');