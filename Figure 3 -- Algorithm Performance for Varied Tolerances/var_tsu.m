% var_tsu.m
% Computes expected posterior covariance using Eq. (3) from Barendregt et
% al., 2022.
function V = var_tsu(T,h0,p)
% Compute integrand:
Int = p.*(exp(-h0*T).*...
    (h0-trapz(h0,h0.*exp(-h0*T).*p)/trapz(h0,exp(-h0*T).*p)).^2+...
    (1-exp(-h0*T)).*...
    (h0-trapz(h0,h0.*(1-exp(-h0*T)).*p)/trapz(h0,(1-exp(-h0*T)).*p)).^2);
% Compute expected posterior covariance:
V = trapz(h0,Int);
end