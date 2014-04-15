% Basic script for visualizing the HDF5 output of the HOS code.

close all


figure(1)
set(gcf,'Position',[50 100 1200 800]);
set(gcf,'Color',[1 1 1]);

SimFolder = '/Users/Claudio/rw/hos/data/demo';


for nfield=[0:1:4000]; 
    
    time = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/time');
    
    eta = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/eta');
    phi = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/phi');
    
    %array1=h5read([RootFolder,'/data_extra',num2str(nfield),'.1.h5'],'/Array1');
    %array2=h5read([RootFolder,'/data_extra',num2str(nfield),'.1.h5'],'/Array2');
    
    Lx = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Lx');
    Ly = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ly');
    Nx = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Nx');
    Ny = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ny');
    
    g = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/g');
   
    Kx = [-Nx/2:1:Nx/2-1]'*2*pi/Lx;
    Ky = [-Ny/2:1:Ny/2-1]'*2*pi/Ly;

    x = Lx*(0:Nx-1)/Nx;
    y = Ly*(0:Ny-1)/Ny;
    
    
    
    figure(1)
    
    col1 = [.7 .9 1];
    col2 = [.7 .9 1];
    ncol1 = 10;
    
    colormap(...
           [ linspace(col1(1),col2(1),ncol1)' linspace(col1(2),col2(2),ncol1)' linspace(col1(3),col2(3),ncol1)' ] ...
    );

    %subplot(2,1,1)
    surf( x, y, eta, ...
          'FaceColor','interp',...
	      'EdgeColor','none',...
	      'FaceLighting','phong'...
        );
    
    
    %axis equal
    daspect([5 5 1])
    axis tight
    camlight left
    view([60,30])
    camzoom(3);
    caxis([-.1 .1]);
    
    set( gca, 'XTick', [], ...
              'YTick', []  ...
       );   
    
   
   figure(2)
   
   colormap('jet');
   Hs = 4*std(reshape(eta,numel(eta),1));
   contour(x, y, eta, [-.5 -.25 .25 .5]*Hs)
   grid on
   
    % Get a frame if you want to create a movie later on.
    % Otherwise comment out the next line.
    Movie(nfield+1) = getframe(gcf);
        
    pause(0.3)
    
end


