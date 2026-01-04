function [Acoef,Adotcoef] = get_coefs2delta_imag_2N(indices,factors,N,kp,delta0_imag)

cscv = @(l,delta0,kp) sin(pi/2*(-l+delta0))./sin(pi/2/N*(-l+delta0)) - sin(pi/2*(2*kp+l+delta0))./sin(pi/2/N*(2*kp+l+delta0));
cscdotv = @(l,delta0,kp) pi/2*cos(pi/2*(-l+delta0))./sin(pi/2/N*(-l+delta0))-pi/2/N*cos(pi/2/N*(-l+delta0)).*sin(pi/2*(-l+delta0))./(sin(pi/2/N*(-l+delta0))).^2+...
    -pi/2*cos(pi/2*(2*kp+l+delta0))./sin(pi/2/N*(2*kp+l+delta0))+pi/2/N*cos(pi/2/N*(2*kp+l+delta0)).*sin(pi/2*(2*kp+l+delta0))./(sin(pi/2/N*(2*kp+l+delta0))).^2;

Acoef = 0; Adotcoef = 0;
for k=1:length(indices)
    Acoef = Acoef + factors(k)*cscv(indices(k),delta0_imag,kp);
    Adotcoef = Adotcoef + factors(k)*cscdotv (indices(k),delta0_imag,kp);
end