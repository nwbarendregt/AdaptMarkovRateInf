% bd_trans_prob.m
% Calculates transition probabilities for a birth-death process using Eq.
% (11) from Barendregt et al., 2022.
function p = bd_trans_prob(mu,lambda,k,i,T,Tp)
% Calculate function parameters:
rho = lambda./mu; a = 2*sqrt(lambda.*mu);
t = T-Tp;
% Set limits of summation:
j = (k+i+2):(k+i+2+50);
% Calculate summation:
s = rho.^(-j/2).*besseli(j,a*t); s(s == Inf) = 0; S = sum(s,'omitnan');
% Calculate transition probabilities:
p = exp(-(lambda+mu)*t).*(rho.^((k-i)/2).*besseli(k-i,a*t)+rho.^((k-i-1)/2).*besseli(k+i+1,a*t));
if S ~= 0
    p = p+exp(-(lambda+mu)*t)*(1-rho).*rho.^k.*S;
end
end