function [ hetaD ] = FtoD( hetaF, hphiF, kx, ky, g )
%   Transforms Fourier coeffs into directional spectrum coeffs.
%
%   [ hetaD ] = FtoD( hetaF, hphiF, kx, ky, g );
%
%   NOTE: the Fourier coeffs are assumed to be
%   ordered with the zero-frequency component in
%   the center of the spectrum. If hu is the output of 
%   fft2, i.e., hu = fft2(u), it is necessary to
%   perform fftshift on hu before using it as input. 


    % Set g=1 unless given as input
    if (nargin==4)
	g=1;
    end
    
    Nx = length(kx);
    Ny = length(ky);
    
    % A loop is highly inefficient, but in this case we prefer to 
    % make the code human readable.
    for ikx = 1:length(kx)

	for iky = 1:length(ky)
	   
	    k = abs( kx(ikx) + 1i*ky(iky) );
	    alpha = sqrt( g/k );
	    
	    hetaD(ikx,iky) = 1i*hphiF(ikx,iky)/alpha + hetaF(ikx,iky);
	     
	    norm = sqrt( g^2/2/(g*k)^(3/2) );
	    hetaF(ikx,iky) = hetaF(ikx,iky)/norm; 
	    hphiF(ikx,iky) = hphiF(ikx,iky)/norm;
	    
	end
    
    end

    % STILL NEED TO NORMALIZE !!!
    
end

