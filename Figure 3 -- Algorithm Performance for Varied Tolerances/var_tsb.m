% var_tsb.m
% Computes expected posterior covariance determinant using Eq. (7) from
% Barendregt et al., 2022.
function V = var_tsb(T,h0,h1,pn,X,Tp)
C = NaN(2,2);
[H0,H1] = meshgrid(h0,h1);
% Calculate transition probabilities conditioned on previous measurement:
if X == 0
    p_0 = H1./(H0+H1)+H0.*exp(-(H0+H1)*(T-Tp))./(H0+H1);
    p_0(1,1) = 1;
    p_1 = H0./(H0+H1)-H0.*exp(-(H0+H1)*(T-Tp))./(H0+H1);
    p_1(1,1) = 0;
else
    p_0 = H1./(H0+H1)-H1.*exp(-(H0+H1)*(T-Tp))./(H0+H1);
    p_0(1,1) = 0;
    p_1 = H0./(H0+H1)+H1.*exp(-(H0+H1)*(T-Tp))./(H0+H1);
    p_1(1,1) = 1;
end
% Compute integrands:
Int_00 = H0-(trapz(h0,H0.*trapz(h1,p_0.*pn.^2,1),2))/...
    trapz(h1,trapz(h0,p_0.*pn.^2,2),1);
Int_01 = H0-(trapz(h0,H0.*trapz(h1,p_1.*pn.^2,1),2))/...
    trapz(h1,trapz(h0,p_1.*pn.^2,2),1);
Int_10 = H1-(trapz(h1,H1.*trapz(h0,p_0.*pn.^2,2),1))/...
    trapz(h1,trapz(h0,p_0.*pn.^2,2),1);
Int_11 = H1-(trapz(h1,H1.*trapz(h0,p_1.*pn.^2,2),1))/...
    trapz(h1,trapz(h0,p_1.*pn.^2,2),1);
% Compute entries of covariance matrix:
C(1,1) = trapz(h1,trapz(h0,Int_00.^2.*p_0.*pn.^2,2),1)+...
    trapz(h1,trapz(h0,Int_01.^2.*p_1.*pn.^2,2),1);
C(2,2) = trapz(h1,trapz(h0,Int_10.^2.*p_0.*pn.^2,2),1)+...
    trapz(h1,trapz(h0,Int_11.^2.*p_1.*pn.^2,2),1);
C(1,2) = trapz(h1,trapz(h0,Int_00.*Int_10.*p_0.*pn.^2,2),1)+...
    trapz(h1,trapz(h0,Int_01.*Int_11.*p_1.*pn.^2,2),1);
C(2,1) = C(1,2);

V = det(C);
end