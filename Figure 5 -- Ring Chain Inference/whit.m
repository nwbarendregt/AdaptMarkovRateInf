% whit.m
% Computes the Whittaker function given by Barendregt et al., 2022.
function W = whit(a,l,m)
% Define integration mesh:
t = linspace(0,10);
[I,J] = size(a); int = NaN(I,J);
for i = 1:I
    for j = 1:J
        % Compute integral component of Whittaker function:
        int(i,j) = trapz(t,t.^(m-l-0.5).*(1+t).^(m+l-0.5).*exp(-a(i,j).*t));
    end
end
% Compute value of Whittaker function:
W = a.^(m+0.5).*exp(-a/2)/gamma(m-l+0.5).*int;