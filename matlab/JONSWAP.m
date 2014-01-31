% Function implementing the JONSWAP directional wave spectrum
%
% [P] = JONSWAP(omega, theta, g, alpha_p, omega_p, gamma, Theta);

function [P] = JONSWAP(omega, theta, g, alpha_p, omega_p, gamma, Theta)

    % default parameters
    if(nargin==2)
	g = 1;
	alpha_p = 2*pi*0.014/10;
        omega_p = 1;
        gamma= 6;
	Theta = 12*pi/180;
    end
		
    
    
    %omega = @(k_x,k_y) sqrt(g*abs(k_x + 1i*k_y));
    %theta = angle(k_x + 1i*k_y);
    
    %[a,b,c]=find(k_x==0 & k_y==0);
    %k_x(a,b)=1;
    %k_y(a,b)=1;
    
    %sigma = @(k_x,k_y) 0.07*( sqrt(g*abs(k_x + 1i*k_y)) < omega_p ) + 0.09*( sqrt(g*abs(k_x + 1i*k_y)) >= omega_p );
    sigma = @(omega) 0.07*( omega < omega_p ) + 0.09*( omega >= omega_p );
    
    P = alpha_p*g^2*omega.^(-5).*exp(-5/4*(omega/omega_p).^(-4)).*...
        gamma.^(exp( -(omega - omega_p).^2))./(2*sigma(omega).^2*omega_p^2).*...
	(2/Theta).*cos(min(abs(pi*theta/Theta),pi/2)).^2;

	
    %P(a,b) = 0;

end