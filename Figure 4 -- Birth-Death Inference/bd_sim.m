% bd_sim.m
% Simulates birth-death process with birth rates lambda and death rates mu.
function [X,T] = bd_sim(Xp,Tp,mu,lambda)
if Xp == 0
    % Calculate transition if X = 0:
    T = Tp-log(rand)/lambda; X = 1;
else
    % Calculate transition if X ~= 0:
    T = Tp-log(rand)/(mu+lambda);
    U = rand;
    if  U < mu/(mu+lambda)
        X = Xp-1;
    else
        X = Xp+1;
    end
end