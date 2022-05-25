% ring_trans_prob.m
% Calculates conditional transition probabilities for ring network using
% Eq. (12) from Barendregt et al., 2022.
function p = ring_trans_prob(A,t,Xp)
% Compute transition probability matrix using matrix exponential:
P = expm(A*t);
% Find conditional transition probability from matrix exponential:
p = P(Xp,:);
end