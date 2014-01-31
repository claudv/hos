close all


figure(1)
set(gcf,'Position',[50 100 1200 800]);
set(gcf,'Color',[1 1 1]);
figure(2)    
set(gcf,'Position',[1400 100 800 1000]);

SimFolder = '/Users/Claudio/rw/hos/data/1';


for nfield=[0:1:4000];
    %sfsafd
    
    time =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/time');
    
    eta =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/eta');
    phi =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/phi');
    
    %array1=h5read([RootFolder,'/data_extra',num2str(nfield),'.1.h5'],'/Array1');
    %array2=h5read([RootFolder,'/data_extra',num2str(nfield),'.1.h5'],'/Array2');
    
    Lx =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Lx');
    Ly =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ly');
    Nx =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Nx');
    Ny =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ny');
    
    Kx = [-Nx/2:1:Nx/2-1]'*2*pi/Lx;
    Ky = [-Ny/2:1:Ny/2-1]'*2*pi/Ly;

    x = Lx*(0:Nx-1)/Nx;
    y = Ly*(0:Ny-1)/Ny;
    
    figure(1)
    
    %subplot(2,1,1)
    surf(x,y,eta,'LineStyle','none')
    xlabel('x','Fontsize',24,'FontName','Times')
    ylabel('y','Fontsize',24,'FontName','Times')
    zlabel('Surface','Fontsize',24,'FontName','Times')
    title(['t=' num2str(time)],'Fontsize',24,'FontName','Times')
    
    set(gca,'Fontsize',24,'FontName','Times')
    
    view([45,60])
    %view([0 0])
    zlim([-1 1]); 
    
%     subplot(2,1,2)
%     mesh(x,y,phi)
%     xlabel('x')
%     ylabel('y')
%     zlabel('phi')
%     title(['t=' num2str(time)])
%     view([45,60])
%     %view([0 0])
%     zlim([-1 1]); 
    
%     subplot(4,1,3)
%     mesh(array1)
%     xlabel('x')
%     ylabel('y')
%     zlabel('array1')
%     title(['t=' num2str(time)])
%     view([45,60])
%     
%     subplot(4,1,4)
%     mesh(array2)
%     xlabel('x')
%     ylabel('y')
%     zlabel('array2')
%     title(['t=' num2str(time)])
%     view([45,60])
    

    Movie(nfield+1) = getframe(gcf);
        
    figure(2)
    
    Spectrum = abs(fftshift(fft2(eta)));
    
    subplot(2,1,1)
    semilogy(Kx,Spectrum(Ny/2+1,:))
    
    subplot(2,1,2)
    imagesc(Kx,Ky,log10(Spectrum))
    %hold on
    %contour(Kx,Ky,log10(Spectrum),[3.8:-0.4:1],'k-')
    caxis([-12 4]);
    axis([-3 3 -3 3]);
    axis equal
    grid on
    colorbar
    
    pause()


   
    
end


