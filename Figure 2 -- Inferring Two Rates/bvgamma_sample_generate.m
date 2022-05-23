% bvgamma_sample_generate.m
% Generates samples from bivariate gamma distribution given in Eq. (8) of
% Barendregt et al., 2022 using the method described in Nadarajah and
% Gupta, 2006.
clear
% Define distribution parameters and number of samples:
m1 = 2; m2 = 2;
a = 1; b = 1; c = a+b;
N_Samples = 1e3;
% Generate bivariate gamma-distributed samples:
W_X = -log(rand(1,N_Samples)); W_Y = -log(rand(1,N_Samples));
W = W_X./(W_X+W_Y);

U = -sum(log(rand(c,N_Samples)),1); U = m1*U;
V = -sum(log(rand(c,N_Samples)),1); V = m2*V;

h0 = U.*W; h1 = V.*W;

save('bvgamma_samples.mat','m1','m2','a','b','h0','h1');