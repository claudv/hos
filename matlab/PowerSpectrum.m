runsubid=1;


T = 50;	% Final simulation time
dtsave = .2;
saveflg = 1;	% 1:standard ouput, 2:extra output

g=1;
delta_k = 0.2;
k_p = 1;
epsilon = 0.1;
sigma = epsilon/sqrt(2)/k_p;

%%% Unscaled power spectral densities
Pk = @(k_x,k_y) sigma^2/(2*pi*delta_k^2)*exp( -0.5*( ( (k_x - k_p)/delta_k/k_p ).^2 + ( k_y/delta_k/k_p ).^2 ) );
%Pk = @(k_x,k_y) JONSWAP(k_x, k_y);


Lx = 128*2*pi/32;
Ly = 128*2*pi/32;
Nx = 4096/32;
Ny = 4096/32;
kx = 2*pi/Lx*[-Nx/2:Nx/2-1];
ky = 2*pi/Ly*[-Ny/2:Ny/2-1];


[KX, KY] = meshgrid(kx,ky);
KA = abs(KX + 1i*KY);

hetar = zeros(Nx,Ny);
%hetar(Nx/2+2,Ny/2+2)=1*Nx*Ny;
rng('shuffle');
hetar = Nx*Ny*sqrt(2*Pk(KX,KY)).*exp(1i*2*pi*rand(Nx,Ny))/16;
%hetar(Ny/2+1+3,Nx/2+1) = 1;
%hetar(Ny/2+1+1,Nx/2+1+4) = 1/2;
hphir = -1i*sqrt(g./KA).*hetar;

hetar(:,1:Nx/2) = 0;
hphir(:,1:Nx/2) = 0;
hetar(1:Ny/2+1,Nx/2+1) = 0;
hphir(1:Ny/2+1,Nx/2+1) = 0;

hetar(Ny/2+1,Nx/2+1) = 0;
hphir(Ny/2+1,Nx/2+1) = 0;

etar = ifft2(ifftshift(hetar));
etar = real(etar);

phir = ifft2(ifftshift(hphir));
phir = real(phir);




rescale = std(reshape(etar,Nx*Ny,1),1)/sigma

etar = etar/rescale;
phir = phir/rescale;

figure()
subplot(1,2,1)
imagesc(etar)
subplot(1,2,2)
imagesc(phir)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write hdf5 datafile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SimFolder = '/Users/Claudio/rw/hos/data/1/';

hdf5File  = [SimFolder, 'initpars.h5'];
unix(['rm ',hdf5File]);

h5create(hdf5File,'/Nx',[1 1]);
h5write(hdf5File, '/Nx', Nx);

h5create(hdf5File,'/Ny',[1 1]);
h5write(hdf5File, '/Ny', Ny);

h5create(hdf5File,'/Lx',[1 1]);
h5write(hdf5File, '/Lx', Lx);

h5create(hdf5File,'/Ly',[1 1]);
h5write(hdf5File, '/Ly', Ly);

h5create(hdf5File,'/T',[1 1]);
h5write(hdf5File, '/T', T);

h5create(hdf5File,'/dtsave',[1 1]);
h5write(hdf5File, '/dtsave', dtsave);

h5create(hdf5File,'/saveflg',[1 1]);
h5write(hdf5File, '/saveflg', saveflg);

h5create(hdf5File,'/runsubid',[1 1]);
h5write(hdf5File, '/runsubid', runsubid);

h5create(hdf5File,'/rampflg',[1 1]);
h5write(hdf5File, '/rampflg', 1);

h5create(hdf5File,'/Tramp',[1 1]);
h5write(hdf5File, '/Tramp', 1);


hdf5File  = [SimFolder, 'initdata.',num2str(runsubid),'.h5'];
unix(['rm ',hdf5File]);

%Matlab is column-major, need to transpose before writing.
h5create(hdf5File,'/eta0',[Ny Nx]);
h5create(hdf5File,'/phi0',[Ny Nx]);

h5write(hdf5File, '/eta0', etar);
h5write(hdf5File, '/phi0', phir);

%h5disp(hdf5File)
