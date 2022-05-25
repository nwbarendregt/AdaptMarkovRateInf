% tsb_sim.m
% Simulates two-state Markov chain with transition rates h0 and h1.
function [X,T] = tsb_sim(Xp,Tp,h0,h1)
if Xp == 0
    % Compute transition time for 0 -> 1 transition:
    dt = -log(rand)/h0;
    X = 1; T = Tp+dt;
else
    % Compute transition time for 1 -> 0 transition:
    dt = -log(rand)/h1;
    X = 0; T = Tp+dt;
end
end