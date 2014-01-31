% The maximum amplidue (amp) for a Stokes wave is around 0.44

function [ eta, phi, phi_x, w, c ] = StokesWave(amp, kap, L, N, t)
     
    g=1;
    dxi = L/N;
    xi = (0:N-1)'*dxi;	

    k = [0:N/2 -N/2+1:-1]'*2*pi/L; 
    sk = sign(k);				
    ka = abs(k);
    
    ome = sqrt(g*kap);
    
    % Use second order solution as the initial guess
    c = ome/kap + 0.5*amp^2*ome^5/kap;
    Stokes_eta = @(x,t)   amp*cos(kap*(x - c*t))  +   0.5*amp^2*abs(kap)*cos(2*kap*(x - c*t));

    eta_guess = Stokes_eta(xi,0);
    
    %options = optimset('Display','iter','TolFun',1e-10);
    %eta = fsolve( @(eta) Fun(eta, c, k, ka), eta_guess, options);

    addpath('./NewtonRaphson 2/');
    options = optimset('TolFun', 1e-12,'TolX',1e-10);
    [eta, resnorm, F, exitflag, output, jacob] = newtonraphson(@(eta) Fun(eta, c, k, ka), eta_guess, options);
    
    %addpath('./NewtonRaphson/');
    %eta = NewtonRaphson(@(eta) Fun(eta, c, k, ka), eta_guess, .1, 100, 'on');
    
    eta = real(ifft(exp(-1i*k*c*t).*fft(eta)));
    heta = fft(eta);
    phi = -c*real(ifft(1i*sk.*heta));
    psi = real(ifft(1i*sk.*fft(phi)));
   
    
    eta_xi = real(ifft(1i*k.*heta));
    phi_xi = real(ifft(1i*k.*fft(phi)));
    psi_xi = real(ifft(1i*k.*fft(psi)));
    x_xi = 1 - real(ifft((-ka.*heta)));

    eta_x = eta_xi./x_xi;
    phi_x = phi_xi./x_xi;
    psi_x = psi_xi./x_xi;
    
    sq = sqrt(1 + eta_x.^2);

    % Velocity component tangent to the surface
    ut = phi_x./sq;
    
    % Velocity component normal to the surface
    un = -psi_x./sq;
    
    % Horizontal and vertical velocity
    u = ut./sq - un.*eta_x./sq;
    w = un./sq + ut.*eta_x./sq;


end
