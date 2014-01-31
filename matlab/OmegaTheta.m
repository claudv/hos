function [ omega, theta ] = OmegaTheta( kx, ky, g )
% Returns polar coordinates on each node of the
% rectangular k-grid in two matrices.
%
% [ omega, theta ] = OmegaTheta( kx, ky, g );
%
% NOTE: since we use meshgrid() we obtain
% size(omega) = [ length(ky) length(kx) ]
% and likewise for theta.

    % Set g=1 unless given as input
    if (nargin==2)
	g=1;
    end
    
    [KX,KY] = meshgrid(kx,ky);
    
    omega = sqrt(g*abs(KX + 1i*KY));
    theta = angle(KX + 1i*KY);
    
    % Enforce 0 < theta <= 2*pi
    theta = theta + 2*pi*(theta < 0);
   

end

