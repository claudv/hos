% Function implementing the JONSWAP directional wave spectrum
%
% [P] = JONSWAP(omega, theta, g, alpha_p, omega_p, gamma, Theta);

function [P] = JONSWAP(omega, theta, g, alpha_p, omega_p, gamma, Theta, theta_shift)

    % default parameters
    if(nargin==2)
	g       = 1;
	alpha_p = 2*pi*0.014/10;
        omega_p = 1;
        gamma   = 6;
	Theta   = 12*pi/180;
	theta_shift = 0;
    end
		
    % convert from 0 < theta < 2*pi to -pi < theta < pi.
    theta = theta - 2*pi*(theta>pi);

    % adding angular shift
    theta = theta + theta_shift;
    
    sigma = @(omega) 0.07*( omega < omega_p ) + 0.09*( omega >= omega_p );
    
    
    P = alpha_p*g^2*omega.^(-5).*exp(-5/4*(omega/omega_p).^(-4)).*...
        gamma.^(exp( -(omega - omega_p).^2./(2*sigma(omega).*sigma(omega)*omega_p^2))).*...
	(2/Theta).*cos(min(abs(pi*theta/Theta),pi/2)).^2;

    
    
    %P(a,b) = 0;
 
end