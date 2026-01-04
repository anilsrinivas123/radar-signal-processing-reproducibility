function [deltahat_real,deltahat_imag] = inverse_mappings(ratio_real,ratio_imag,N,max_frequency_bins,delta_ini,maxiter)

delta0_real = delta_ini; delta0_imag = delta_ini; 
for iter=1:maxiter
    [A1,A1dot] = get_coefs2delta_real_2N([-1 0 1],[-1 0 1],N,max_frequency_bins-1,delta0_real);
    [B1,B1dot] = get_coefs2delta_real_2N([-1 0 1],[-1 2 -1],N,max_frequency_bins-1,delta0_real);
    const_real = A1./B1;
    cN_real = (A1dot.*B1 - A1.*B1dot)./B1.^2;
    deltahat_real = (ratio_real-const_real)./cN_real + delta0_real;
    %est_real(iter,:) = deltahat_real;
    delta0_real = deltahat_real;
    
    [A2,A2dot] = get_coefs2delta_imag_2N([-1 0 1],[-1 0 1],N,max_frequency_bins-1,delta0_imag);
    [B2,B2dot] = get_coefs2delta_imag_2N([-1 0 1],[-1 2 -1],N,max_frequency_bins-1,delta0_imag);
    const_imag= A2./B2;
    cN_imag= (A2dot.*B2 - A2.*B2dot)./B2.^2;
    deltahat_imag = (ratio_imag-const_imag)./cN_imag + delta0_imag;
    %est_imag(iter,:) = deltahat_imag;
    delta0_imag = deltahat_imag;
end
end