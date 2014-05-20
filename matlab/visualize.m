% Basic script for visualizing the HDF5 output of the HOS code.

close all



SimFolder = '/Users/Claudio/rw/hos/data/demo';


for nfield=[0:1:4000];
    
    
    time =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/time');
    
    eta =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/eta');
    phi =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/phi');
    
    %array1=h5read([SimFolder,'/data_extra',num2str(nfield),'.1.h5'],'/Array1');
    array2=h5read([SimFolder,'/data_extra',num2str(nfield),'.1.h5'],'/Array2');
    
    Lx =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Lx');
    Ly =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ly');
    Nx =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Nx');
    Ny =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ny');
    
    g = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/g');
   
    Kx = [-Nx/2:1:Nx/2-1]'*2*pi/Lx;
    Ky = [-Ny/2:1:Ny/2-1]'*2*pi/Ly;

    x = Lx*(0:Nx-1)/Nx;
    y = Ly*(0:Ny-1)/Ny;
    
    std(reshape(eta,Nx*Ny,1),1)
    
    
    figure(1)
    set(gcf,'Position',[50 100 1200 800]);
    set(gcf,'Color',[1 1 1]);
  

    %subplot(2,1,1)
    surf(x,y,eta,array2,'FaceColor','interp',...
	         'EdgeColor','none',...
	         'FaceLighting','phong')
    xlabel('x','Fontsize',24,'FontName','Times')
    ylabel('y','Fontsize',24,'FontName','Times')
    zlabel('Surface','Fontsize',24,'FontName','Times')
    title(['t=' num2str(time)],'Fontsize',24,'FontName','Times')
    
    set(gca,'Fontsize',24,'FontName','Times')
    %axis equal
    daspect([5 5 1])
    axis tight
    camlight left
    view([25,40])
    %caxis([-.1 .1]);
    
    
    % Get a frame if you want to create a movie later on.
    % Otherwise comment out the next line.
    %Movie(nfield+1) = getframe(gcf);
        
    figure(2)    
    set(gcf,'Position',[1400 100 800 1000]);
    
   
    [ hetaD ] = FtoD( fftshift(fft2(eta.')), fftshift(fft2(phi.')), Kx, Ky, g );
    SpectrumD = abs(hetaD);
    SpectrumD = SpectrumD.';
    SpectrumD = SpectrumD/max(max(abs(SpectrumD)));
   
    Spectrum = abs(fftshift(fft2(eta)));
    %Spectrum0 = Spectrum;
    Spectrum = Spectrum/max(max(abs(Spectrum)));
    
    subplot(3,1,1)
    loglog(Kx,Spectrum(Ny/2+1,:))
    
    subplot(3,1,2)
    imagesc(Kx,Ky,log10(Spectrum))
    set(gca,'YDir','normal');
    %hold on
    %contour(Kx,Ky,log10(Spectrum),[3.8:-0.4:1],'k-')
    caxis([-3 0]);
    axis([-12 12 -12 12]);
    axis equal
    grid on
    colorbar
    
    subplot(3,1,3)
    [ omega, theta ] = OmegaTheta( Kx, Ky, g);
    [X,Y] = pol2cart(theta,omega);
    contourf(X,Y,SpectrumD,[.1 .2 .3 .4 .5 .6 .7 .8 .9],'EdgeColor','None');
    axis([-10 10 -10 10]);
    %axis equal
    axis square
    grid on
    colorbar
    
    pause()


   
    
end


