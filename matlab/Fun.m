function [ f ] = Fun(eta, c, k, ka)
    
    heta = fft(eta);
    
    heta_xi = 1i*k.*heta;
    
    eta_xi = real(ifft(heta_xi));
    chi_xi = 1 + real(ifft(ka.*heta));
  
    J = chi_xi.^2 + eta_xi.^2;
    
    f = 0.5*c^2.*( 1./J - 1) + eta;

end

