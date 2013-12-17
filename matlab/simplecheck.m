close all


figure(1)
set(gcf,'Position',[50 100 600 700]);
figure(2)    
set(gcf,'Position',[800 100 600 700]);

DataFolder1 = '/Users/Claudio/rw/hos/data/1';
DataFolder2 = '/Users/Claudio/rw/hos/data/2';
SubID1 = '1';
SubID2 = '1';

for nfield=[0:1:4000];
    
    
    time =h5read([DataFolder1,'/data',num2str(nfield),'.',SubID1,'.h5'],'/time');
    time2=h5read([DataFolder2,'/data',num2str(nfield),'.',SubID2,'.h5'],'/time');
    
    eta =h5read([DataFolder1,'/data',num2str(nfield),'.',SubID1,'.h5'],'/eta');
    eta2=h5read([DataFolder2,'/data',num2str(nfield),'.',SubID2,'.h5'],'/eta');
    phi =h5read([DataFolder1,'/data',num2str(nfield),'.',SubID1,'.h5'],'/phi');
    phi2=h5read([DataFolder2,'/data',num2str(nfield),'.',SubID2,'.h5'],'/phi');
    
    %array1_1=h5read([DataFolder1,'/data_extra',num2str(nfield),'.1.h5'],'/Array1');
    %array1_2=h5read([DataFolder2,'/data_extra',num2str(nfield),'.1.h5'],'/Array1');
    
    Lx =h5read([DataFolder1,'/data',num2str(nfield),'.',SubID1,'.h5'],'/Lx');
    Ly =h5read([DataFolder1,'/data',num2str(nfield),'.',SubID1,'.h5'],'/Ly');
    Nx =h5read([DataFolder1,'/data',num2str(nfield),'.',SubID1,'.h5'],'/Nx');
    Ny =h5read([DataFolder1,'/data',num2str(nfield),'.',SubID1,'.h5'],'/Ny');
    
    Kx = [-Nx/2+1:1:Nx/2]'*2*pi/Lx;
    Ky = [-Ny/2+1:1:Ny/2]'*2*pi/Ly;

    err(nfield+1) = norm(eta-eta2);
    
    figure(1)
    
    subplot(4,1,1)
    mesh(eta)
    xlabel('x')
    ylabel('y')
    zlabel('eta')
    title(['t=' num2str(time)])
    view([45,60])
    %view([0 0])
    %zlim([-.2 .2]); 
    
    subplot(4,1,2)
    mesh(eta2)
    xlabel('x')
    ylabel('y')
    zlabel('phi')
    title(['t=' num2str(time)])
    view([45,60])
    %view([0 0])
    %zlim([-.2 .2]); 
    
    subplot(4,1,3)
    mesh(phi)
    xlabel('x')
    ylabel('y')
    zlabel('array1')
    title(['t=' num2str(time2)])
    view([45,60])
    
    subplot(4,1,4)
    mesh(phi2)
    xlabel('x')
    ylabel('y')
    zlabel('array2')
    title(['t=' num2str(time2)])
    view([45,60])
    
    %set(gcf,'Position',[800 800 600 600]);
    
    
    figure(2)
    
    subplot(2,1,1)
    semilogy(abs(fft2(eta(1,:))))
    
    subplot(2,1,2)
    %mesh(abs(fft2(eta)))
    %set(gca,'ZScale','log')
    %contour(log(abs(fftshift(fft2(eta-etat.')))))
    imagesc(Kx,Ky,log(abs(fftshift(fft2(eta)))))
    %caxis([-9 0]);
    %axis([-5 5 -2 2]);
    axis equal
    grid on
    colorbar
    
    pause()

%     for ix=1:Nx;
% 	
% 	for iy=1:Ny
% 
% 		x = dxi*(ix-1);
% 		y = dyi*(iy-1);
% 		
% 		etat(ix,iy) = Stokes_2d_eta(x,y,time);
% 		phit(ix,iy) = Stokes_2d_phi(x,y,time);
% 	
% 	end
% 	
%     end
%     
%     subplot(3,1,3)
%     mesh(phit.')
%     xlabel('x')
%     ylabel('y')
%     zlabel('phi')
%     title(['t=' num2str(time)])
%     view([45,60])
%     zlim([-.2 .2]); 
%     
%     errn(nfield+1) = max(max(abs(phi.' - phit, Inf)))/amp;
   
    
end


