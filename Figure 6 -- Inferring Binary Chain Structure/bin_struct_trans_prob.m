% bin_struct_trans_prob.m
% Calculates conditional transition probabilities using matrix exponential
% of infinitesimal generator matrix.
function p = bin_struct_trans_prob(A,t,X,Xp)
% Calculate general transition probability matrix:
p = expm(A*t);
% Calculate conditional transition probability from general matrix:
p = p(Xp,:);
p = p(X);
end