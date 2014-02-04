runid=1;
runsubid=1;

%Lx=32*2*pi;
%Ly=32*2*pi;
Lx=4*pi;
Ly=4*pi/sqrt(3);

Nx=12*64;
Ny=4*64;
%Nx=2*32*64;
%Ny=2*32*64;
N=Nx*Ny;

dx = Lx/Nx;
x = (0:Nx-1)'*dx;

%kx = [0:Nx/2 - Nx/2+1:-1]'*pi/Lx; 

dy = Ly/Ny;
y = (0:Ny-1)'*dy;	

%ky = [0:Ny/2 - Ny/2+1:-1]'*pi/Ly; 

% Final simulation time
T=20;
dtsave = 1;

amp = 0.15;

%kx = 32*2*pi/Lx;
%ky = 0*2*pi/Ly;
kx = 1/2;
ky = sqrt(3)/2;
%kx = 2*pi/Lx;
%ky = 0;


kap = sqrt(kx^2+ky^2);

ome = sqrt(kap);

% c = ome/kap + 0.5*amp^2*ome^5/kap;	
% Stokes_2d_eta = @(x,y,t)   amp*cos(kx*x + ky*y - kap*c*t)   + 0.5*amp^2*abs(kap)*cos(2*kx*x + 2*ky*y - 2*kap*c*t);
% Stokes_2d_phi = @(x,y,t) c*amp*sin(kx*x + ky*y - kap*c*t) + 0.5*c*amp^2*abs(kap)*sin(2*kx*x + 2*ky*y - 2*kap*c*t);


N = 2048;
L = ceil( sqrt( Lx^2 + Ly^2 )/(2*pi/kap) )*2*pi/kap;

[ etas, phis, phis_xs, ws, c_stokes ] = StokesWave(amp, kap, L, N, 0);

hetas = fft(etas);
hphis = fft(phis);
hws = fft(ws);

kappa = [0:N/2 -N/2+1:-1]'*2*pi/L; 
sk = sign(kappa);
xis = (0:N-1).'*L/N;
xs = xis - real(ifft(( 1i*sk.*hetas )));


Nlarge = 32*N;
xis_large = (0:Nlarge-1).'*L/Nlarge;

hetas_pad = zeros(Nlarge,1);
hetas_pad(1:N/2+1) = hetas(1:N/2+1)*Nlarge/N; 
hetas_pad(end-N/2+2:end) = hetas(N/2+2:end)*Nlarge/N;
etas_int = real(ifft(hetas_pad));

hphis_pad = zeros(Nlarge,1);
hphis_pad(1:N/2+1) = hphis(1:N/2+1)*Nlarge/N; 
hphis_pad(end-N/2+2:end) = hphis(N/2+2:end)*Nlarge/N;
phis_int = real(ifft(hphis_pad));

hws_pad = zeros(Nlarge,1);
hws_pad(1:N/2+1) = hws(1:N/2+1)*Nlarge/N; 
hws_pad(end-N/2+2:end) = hws(N/2+2:end)*Nlarge/N;
ws_int = real(ifft(hws_pad));

hxs_pad = zeros(Nlarge,1);
hxs = fft(xs - xis);
hxs_pad(1:N/2+1) = hxs(1:N/2+1)*Nlarge/N; 
hxs_pad(end-N/2+2:end) = hxs(N/2+2:end)*Nlarge/N;
xs_int = real(ifft(hxs_pad)) + xis_large;


Stokes_eta_exact = @(xx) interp1([ xs_int  - L;  xs_int; xs_int  + L], [etas_int; etas_int; etas_int], xx, 'cubic');
Stokes_phi_exact = @(xx) interp1(xs_int, phis_int, xx, 'cubic');
Stokes_w_exact = @(xx) interp1(xs_int, ws_int, xx, 'cubic');

etasx = Stokes_eta_exact(xis);
phisx = Stokes_phi_exact(xis);
wsx = Stokes_w_exact(xis);

Stokes_2d_eta = @(x,y,t) Stokes_eta_exact( (kx*x + ky*y)/kap );
Stokes_2d_phi = @(x,y,t) Stokes_phi_exact( (kx*x + ky*y)/kap );
Stokes_2d_w   = @(x,y,t) Stokes_w_exact( (kx*x + ky*y)/kap );

eta0 = zeros(Nx,Ny);
phi0 = zeros(Nx,Ny);
w0 = zeros(Nx,Ny);

%for ix=1:Nx;
    for ny=1:Ny

        %x = dxi*(ix-1)
        y = dy*(ny-1);

        eta0(:,ny) = Stokes_2d_eta(x,y,0) + 0.1*cos(kx*x); %...
			%+ 0.001*cos(x/8) + 0.001*cos(x/4) + 0.001*cos(x/2) ...
			%+ 0.001*cos(x/8 + y/8) + 0.001*cos(x/4 + y/4) + 0.001*cos(x/4 + y/2);
        phi0(:,ny) = Stokes_2d_phi(x,y,0);
        w0(:,ny) = Stokes_2d_w(x,y,0);
        %eta0(ix,iy) = 0.1*sin(kx*dxi*(ix-1));
        %phi0(ix,iy) = 0.05/sqrt(sqrt(kx^2+ky^2))*cos(kx*dxi*(ix-1) + ky*dyi*(iy-1));
        %phi0(ix,iy) = cos(2*pi*dyi*(iy-1)/Ly);
        %phi0(ix,iy) = ix + Nx*(iy-1);
    
    end
%end

% Add perturbation
%hpert = 0.005*rand(size(eta0)) + 0.005*1i*rand(size(eta0));
%hpert(1,1) = 0;
%hpert = hpert-mean(mean(hpert));
%eta0 = eta0 + real(ifft2(hpert));




hdf5File  = ['../',num2str(runid),'/initpars.h5'];
unix(['rm ',hdf5File]);

h5create(hdf5File,'/Nx',[1 1]);
h5write(hdf5File, '/Nx', Nx);

h5create(hdf5File,'/Ny',[1 1]);
h5write(hdf5File, '/Ny', Ny);

h5create(hdf5File,'/Lx',[1 1]);
h5write(hdf5File, '/Lx', Lx);

h5create(hdf5File,'/Ly',[1 1]);
h5write(hdf5File, '/Ly', Ly);

h5create(hdf5File,'/g',[1 1]);
h5write(hdf5File, '/g', g);

h5create(hdf5File,'/T',[1 1]);
h5write(hdf5File, '/T', T);

h5create(hdf5File,'/dtsave',[1 1]);
h5write(hdf5File, '/dtsave', dtsave);

h5create(hdf5File,'/saveflg',[1 1]);
h5write(hdf5File, '/saveflg', 2);

h5create(hdf5File,'/runsubid',[1 1]);
h5write(hdf5File, '/runsubid', runsubid);

h5create(hdf5File,'/rampflg',[1 1]);
h5write(hdf5File, '/rampflg', 1);

h5create(hdf5File,'/Tramp',[1 1]);
h5write(hdf5File, '/Tramp', 60);



hdf5File  = [num2str(runid),'/initdata.',num2str(runsubid),'.h5'];
unix(['rm ',hdf5File]);

%Matlab is column-major, need to transpose before writing.
h5create(hdf5File,'/eta0',[Ny Nx]);
h5create(hdf5File,'/phi0',[Ny Nx]);

h5write(hdf5File, '/eta0', eta0.');
h5write(hdf5File, '/phi0', phi0.');

%h5disp(hdf5File)