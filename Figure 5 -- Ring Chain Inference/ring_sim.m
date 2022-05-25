% ring_sim.m
% Simulates ring chain of states with clockwise transition rate h_+ and
% counter-clockwise transition rate h_-.
function [X,T] = ring_sim(Xp,Tp,hp,hm,m)
    % Draw transition time:
    T = Tp-log(rand)/(hp+hm);
    % Draw transition direction:
    U = rand;
    if  U < hp/(hp+hm)
        X = (Xp+1)*(Xp < m)+1*(Xp == m);
    else
        X = (Xp-1)*(Xp > 1)+m*(Xp == 1);
    end
end