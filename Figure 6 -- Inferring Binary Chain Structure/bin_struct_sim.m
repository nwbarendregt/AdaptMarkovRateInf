%bin_struct_sim.m
% Simulates arbitrary Markov chain given infinitesimal generator matrix A:
function [X,T] = bin_struct_sim(Xp,Tp,A)
% Draw transition time:
T = Tp+log(rand)/A(Xp,Xp);
% Find possible transitions:
a = A(Xp,:); a = a(1:length(a) ~= Xp);
% Draw transition direction:
U = rand; X = 1; Phi = -a(X)/A(Xp,Xp);
    while Phi < U
        X = X+1;
        Phi = Phi-a(X)/A(Xp,Xp);
    end
end