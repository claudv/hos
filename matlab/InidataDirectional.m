% Simple script that generates HDF5 input data
% files for HOS simulations. It creates the files:
%
% * initpars.h5       : simulation parameters
% * initdata.runid.h5 : initial condition (eta and phi fields)
%
% The initial condition is created from a
% superposition of Fourier modes with prescribed
% directional specrtum and random phase.
clear all
close all


T = 150;			% Final simulation time
dtsave =  .04;		% Save snapshots every dtsave interval
saveflg = 2;		% 1:standard ouput, 2:extra output
g=9.8;			% Gravity constant

runsubid=1;		% You can label your run with this ID

Lx = 4*2*pi;	% Domain size in x
Ly = 4*2*pi;	% Domain size in y

Nx = 2*512;		% Number of Fourier modes in x
Ny = 2*512;		% Number of Fourier modes in y

kx = 2*pi/Lx*[-Nx/2:Nx/2-1];
ky = 2*pi/Ly*[-Ny/2:Ny/2-1];

    


% Setting free surface variance.
Hs = 0.08;
sigma = Hs/4;


[ omega, theta ] = OmegaTheta( kx, ky, g );

% JONSWAP spectrum UNSCALED!
% First spectral peak
alpha_p = 2*pi*0.014/10;
omega_p = 2*pi;
gamma = 6;
Theta = 12*pi/180;
	
Pdir = JONSWAP(omega, theta, g, alpha_p, omega_p, gamma, Theta, 0);

% Second spectral peak
% alpha_p = 2*pi*0.014/10;
% omega_p = 2*pi;
% gamma = 6;
% Theta = 45*pi/180;
% 	
% Pdir = Pdir + 0.5*JONSWAP(omega, theta, g, alpha_p, omega_p, gamma, Theta, pi/4);

rng('shuffle');
hetaD = sqrt(Pdir.').*exp(1i*2*pi*rand(Nx,Ny));


[ hetar, hphir ] = DtoF( hetaD, kx, ky, g );
cutoff = ones(size(omega)).*(omega < 3*omega_p);    
hetar = cutoff.*hetar;
hphir = cutoff.*hphir;


% PROBABLY NEED TO TRANSPOSE AT THE END
hetar = hetar.';
hphir = hphir.';


figure
imagesc(kx,ky,log10(abs(hetar)))
%caxis(max(max(log10(abs(hetar))))*[ .01 0.9 ]);
set(gca,'YDir','normal');


figure
[X,Y] = pol2cart(theta,omega);
%contourf(X,Y,Pdir,[cont1:step:maxs],'EdgeColor','None')
mesh(X,Y,Pdir)
	
% Compute wave field in the physical space
etar = ifft2(ifftshift(hetar));
%etar = real(etar);

phir = ifft2(ifftshift(hphir));
%phir = real(phir);

% Rescale the wave field to the desired eta variance
rescale = std(reshape(etar,Nx*Ny,1),1)/sigma;
etar = etar/rescale;
phir = phir/rescale;

% Plot the wave field
figure()
subplot(1,2,1)
imagesc([0 Lx],[0 Ly], etar)
set(gca,'YDir','normal');
subplot(1,2,2)
imagesc([0 Lx],[0 Ly], phir)
set(gca,'YDir','normal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write hdf5 datafile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SimFolder = '/Users/Claudio/rw/hos/data/demo/';
%SimFolder = '/Users/Claudio/rw/hos/Csource/2dpar/';

hdf5File  = [SimFolder, 'initpars.h5'];
unix(['rm ', hdf5File]);

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
h5write(hdf5File, '/saveflg', saveflg);

h5create(hdf5File,'/runsubid',[1 1]);
h5write(hdf5File, '/runsubid', runsubid);

h5create(hdf5File,'/rampflg',[1 1]);
h5write(hdf5File, '/rampflg', 0);

h5create(hdf5File,'/Tramp',[1 1]);
h5write(hdf5File, '/Tramp', 60);

h5create(hdf5File,'/windflg',[1 1]);
h5write(hdf5File, '/windflg', 0);

hdf5File  = [SimFolder, 'initdata.',num2str(runsubid),'.h5'];
unix(['rm ',hdf5File]);

%Matlab is column-major, need to transpose before writing.
h5create(hdf5File,'/eta0',[Ny Nx]);
h5create(hdf5File,'/phi0',[Ny Nx]);

h5write(hdf5File, '/eta0', etar);
h5write(hdf5File, '/phi0', phir);

%h5disp(hdf5File)
