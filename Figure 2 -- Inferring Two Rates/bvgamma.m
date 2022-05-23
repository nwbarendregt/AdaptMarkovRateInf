% bvgamma.m
% Defines bivariate gamma distribution function given in Eq. (8) of
% Barendregt et al., 2022.
function f = bvgamma(x,y,a,b,m1,m2)
% Compute model parameters:
c = a+b; C = ((m1*m2)^c*gamma(c)*gamma(a)*gamma(b))^(-1);
% Compute value of distribution:
f = C*gamma(b)*(x.*y).^(c-1).*(x/m1+y/m2).^((a-1)/2-c).*...
    exp(-0.5*(x/m1+y/m2)).*whit(x/m1+y/m2,c-a+(1-a)/2,c-a/2);
end