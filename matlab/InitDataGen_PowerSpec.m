% Simple script that generates HDF5 input data
% files for HOS simulations. It creates the files:
%
% * initpars.h5       : simulation parameters
% * initdata.runid.h5 : initial condition (eta and phi fields)
%
% The initial condition is created from a
% superposition of Fourier modes with prescribed
% amplitude and random phase.
clear all
close all


T = 10;			% Final simulation time
dtsave = .04;		% Save snapshots every dtsave interval
saveflg = 2;		% 1:standard ouput, 2:extra output
g=9.8;			% Gravity constant

runsubid=1;		% You can label your run with this ID

Lx = 4*pi;	% Domain size in x
Ly = 4*pi;	% Domain size in y

Nx = 512;		% Number of Fourier modes in x
Ny = 512;		% Number of Fourier modes in y

kx = 2*pi/Lx*[-Nx/2:Nx/2-1];
ky = 2*pi/Ly*[-Ny/2:Ny/2-1];

SpectrumType = 'Fourier';
%SpectrumType = 'Directional';

switch SpectrumType
    
    case('Fourier')
    
	% Gaussian spectrum UNSCALED!
	k_p = 1;
	delta_k = 0.01;
	epsilon = 0.1;
	sigma = epsilon/sqrt(2)/k_p;

	Pk = @(k_x,k_y) sigma^2/(2*pi*delta_k^2)*exp( -0.5*( ( (k_x - k_p)/delta_k/k_p ).^2 + ( k_y/delta_k/k_p ).^2 ) );

	[KX, KY] = meshgrid(kx,ky);
	KA = abs(KX + 1i*KY);

    % good practice to reshuffle rng to avoid bad suprises
	rng('shuffle');


	% SURFACE ELEVATION -----------------------------%
	hetar = sqrt(Pk(KX,KY)).*exp(1i*2*pi*rand(Ny,Nx));
    %hetar = zeros(Ny,Nx);
    %hetar(Ny/2+1+2,Nx/2+1+1)=1;
    %hetar(Ny/2+1+2,Nx/2+1+2)=0.0001;
    
	% put Nyquist freqs equal zero
	hetar(:,1) = 0;
	hetar(1,:) = 0;

	% enforce zero mean
	hetar(Ny/2+1,Nx/2+1) = 0;

	% make it Hermitian
	hetar(2:Ny/2,Nx/2+1) = conj(flipud(hetar(Ny/2+2:end,Nx/2+1)));
	hetar(2:end,2:Nx/2)  = conj(rot90(hetar(2:end,Nx/2+2:end),2));
	%------------------------------------------------%

	% VELOCITY POTENTIAL ----------------------------%
	% We choose right-going waves
	hphir = zeros(Ny,Nx);

	% fill in half elements
	hphir(2:end,Nx/2+2:end)  = -1i*sqrt(g./KA(2:end,Nx/2+2:end)).*hetar(2:end,Nx/2+2:end);
	hphir(Ny/2+2:end,Nx/2+1) = -1i*sqrt(g./KA(Ny/2+2:end,Nx/2+1)).*hetar(Ny/2+2:end,Nx/2+1);

	% make it Hermitian
	hphir(2:Ny/2,Nx/2+1) = conj(flipud(hphir(Ny/2+2:end,Nx/2+1)));
	hphir(2:end,2:Nx/2)  = conj(rot90(hphir(2:end,Nx/2+2:end),2));

	hphir(Ny/2+1,Nx/2+1) = 0;
    %------------------------------------------------%

    case('Directional')

	% Setting free surface variance.
	sigma = 8.38/4;
	
	% JONSWAP spectrum UNSCALED!
	%[ omega, theta ] = OmegaTheta( kx, ky );
	%Pk = JONSWAP(omega, theta);
	%rng('shuffle');
	%hetaD = sqrt(Pk.').*exp(1i*2*pi*rand(Nx,Ny));
	
	% 
	%hetaD = zeros(Nx,Ny);
	%hetaD(66,33)=1;
	%hetaD(64,33)=0;
	
	% Compute polar coordinates on the rectangular Fourier
	% grid.
	[ omega, theta ] = OmegaTheta( kx, ky, g );
	omegamax = max(max(omega));
	
	% Read directional spectrum data from spec
	% file.
	addpath('/Users/Claudio/rw/hos/matlab/andrea');
	date='9 November 2007 00:00:00';
	numdate=datenum(date);
	specdata=fun_2Dspec_read('andrea/andrea.spec',1,1);
	specall=specdata.SPEC;
	datevec=specdata.DATE;
	itime=find(datevec==numdate);
	ntime=length(datevec)
	
	

	freq=specdata.FREQ;
	dir=specdata.DIR;
	nfreq=length(freq)
	
	% Extract the spectrum corresponding to
	% the chosen time.
	spec=specall((itime-1)*nfreq+1:itime*nfreq,:);
	
	% Filtering high frequencies
	spec(end-8:end,:) = 0;
	
	
	% Convert angles from oceanographic convention 
	% into rest-of-the-world convention.
	dir = pi/2 - dir;
	dir = dir + 2*pi*(dir < 0);
	
	% Sort angles in ascending order
	[dir,dir_idx]=sort(dir');
	dir = dir';
	spec_sort=zeros(size(spec));
	for i=1:length(dir)
		spec_sort(:,i)=spec(:,dir_idx(i));
	end
	spec = spec_sort;
	clear spec_sort
	
	
	% Padding spec data to avoid troubles
	% during interpolation, by adding points at freq = 0
	% and theta = 2*pi.
	% Also, padding beyond freq > omegamax/2/pi.
	dfreq = freq(end)-freq(end-1);
	freqpad = ( freq(end) + dfreq : dfreq : omegamax/2/pi + dfreq )';
	freq = [ 0; freq; freqpad ];
	dir  = [ dir(end)-2*pi; dir; dir(1) + 2*pi ];
	spec = [ spec(:,end) spec spec(:,1) ];
	spec = [ zeros(1,length(dir)) ; spec; zeros(length(freqpad),length(dir))];
	
	
	
	[ Omega, Theta ] = meshgrid(2*pi*freq, dir);
	
	figure
	cont1=0.1;maxs=100;
	step=(maxs-cont1)/20;
	[X,Y] = pol2cart(Theta,Omega);
	%contourf(X,Y,spec',[cont1:step:maxs],'EdgeColor','None')
	mesh(X,Y,spec')
	%spec = spec.';
	
	% Interpolate spectrum on the rectangular
	% grid.
	Pdir = griddata(Omega,Theta,spec',omega,theta,'cubic');
	
	figure
	[X,Y] = pol2cart(theta,omega);
	%contourf(X,Y,Pdir,[cont1:step:maxs],'EdgeColor','None')
	mesh(X,Y,Pdir)
	
	hetaD = sqrt(Pdir.').*exp(1i*2*pi*rand(Nx,Ny));
	
	
	[ hetar, hphir ] = DtoF( hetaD, kx, ky, g );
	
	% PROBABLY NEED TO TRANSPOSE AT THE END
	hetar = hetar.';
	hphir = hphir.';

end

figure
imagesc(kx,ky,log10(abs(hetar)))
set(gca,'YDir','normal');
    
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
h5write(hdf5File, '/windflg', 1);

h5create(hdf5File,'/Uwind_x',[1 1]);
h5write(hdf5File, '/Uwind_x', 1.1);

h5create(hdf5File,'/Uwind_y',[1 1]);
h5write(hdf5File, '/Uwind_y', 0.4);


hdf5File  = [SimFolder, 'initdata.',num2str(runsubid),'.h5'];
unix(['rm ',hdf5File]);

%Matlab is column-major, need to transpose before writing.
h5create(hdf5File,'/eta0',[Ny Nx]);
h5create(hdf5File,'/phi0',[Ny Nx]);

h5write(hdf5File, '/eta0', etar);
h5write(hdf5File, '/phi0', phir);

%h5disp(hdf5File)
