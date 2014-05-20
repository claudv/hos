close all



SimFolder = '/Users/Claudio/rw/hos/data/demo';

nfield=0;
    
    
time =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/time');

eta =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/eta');
phi =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/phi');

eta_t=h5read([SimFolder,'/data_extra',num2str(nfield),'.1.h5'],'/Array1');

Lx =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Lx');
Ly =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ly');
Nx =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Nx');
Ny =h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ny');

g = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/g');

Kx = [-Nx/2:1:Nx/2-1]'*2*pi/Lx;
Ky = [-Ny/2:1:Ny/2-1]'*2*pi/Ly;
[KX, KY] = meshgrid(Kx,Ky);

heta = fftshift(fft2(eta));
heta_x = 1i*KX.*heta;
heta_x = ifftshift(heta_x);
heta = ifftshift(heta);
eta_x = real(ifft2(heta_x));

heta_t = fft2(eta_t);

%DeltaK = 10;
%hW=exp(-(KX.^2 + KY.^2)/DeltaK^2);
%hW=ifftshift(hW);
%W = real(ifft2(hW));
%W=ones(size(eta));

% eta aurocorrelation function
W=real(ifft2(heta.*conj(heta)));

%W=real(ifft2(heta_x.*conj(heta_t)));

% eta smoothed ----------------------------------------%
%KX_shift = ifftshift(KX);
%KY_shift = ifftshift(KY);
%eta_smooth = real(ifft2(heta.*((KX_shift.^2+KY_shift.^2)<100)));
%eta_t_smooth = real(ifft2(heta_t.*((KX_shift.^2+KY_shift.^2)<100)));
hWW = fft2(W.^2);
eta_smooth = real(ifft2(heta.*conj(hWW)))/hWW(1,1);
eta_t_smooth = real(ifft2(heta_t.*conj(hWW)))/hWW(1,1);

[ CPlus, CMinus, LambdaPlus, LambdaMinus, XPlusx, XPlusy, XMinusx, XMinusy ] = LocalPhaseSpeed( eta, eta_t, W, Kx, Ky );
Cx = CPlus.*XPlusx;
Cy = CPlus.*XPlusy;

x = Lx*(0:Nx-1)/Nx;
y = Ly*(0:Ny-1)/Ny;

std(reshape(eta,Nx*Ny,1),1)

figure(1)
set(gcf,'Position',[50 100 1200 800]);
set(gcf,'Color',[1 1 1]);

surf(x,y,eta,'FaceColor','interp',...
         'EdgeColor','none',...
         'FaceLighting','phong')
xlabel('x','Fontsize',24,'FontName','Times')
ylabel('y','Fontsize',24,'FontName','Times')
zlabel('Surface','Fontsize',24,'FontName','Times')
title(['t=' num2str(time)],'Fontsize',24,'FontName','Times')

colormap('Gray');
set(gca,'Fontsize',24,'FontName','Times')
%axis equal
daspect([5 5 1])
axis tight
camlight left
view([15,40])
caxis([-.1 .1]);



figure(2)
set(gcf,'Position',[1400 100 800 1000]);

subplot(2,1,1)
surf(x,y,Cx,'EdgeColor','none')
xlabel('x','Fontsize',24,'FontName','Times')
ylabel('y','Fontsize',24,'FontName','Times')
zlabel('Cx','Fontsize',24,'FontName','Times')
title(['t=' num2str(time)],'Fontsize',24,'FontName','Times')

set(gca,'Fontsize',24,'FontName','Times')
axis tight
view([15,40])

subplot(2,1,2)
surf(x,y,Cy,'EdgeColor','none')
xlabel('x','Fontsize',24,'FontName','Times')
ylabel('y','Fontsize',24,'FontName','Times')
zlabel('Cy','Fontsize',24,'FontName','Times')
title(['t=' num2str(time)],'Fontsize',24,'FontName','Times')

set(gca,'Fontsize',24,'FontName','Times')
axis tight
view([15,40])


figure(3)

%set(cga,'FontSize',15);

stride=40;
[X, Y] = meshgrid(x(1:stride:end),y(1:stride:end));

subplot(2,2,1)
set(gca,'FontSize',16)
imagesc([0 Lx],[0 Ly], CPlus)
colorbar
set(gca,'YDir','normal');
hold on
quiver(X,Y,XPlusx(1:stride:end,1:stride:end),...
           XPlusy(1:stride:end,1:stride:end),1)
title('c^+');
       
       
subplot(2,2,2)
set(gca,'FontSize',16)
imagesc([0 Lx],[0 Ly], CMinus)
colorbar
set(gca,'YDir','normal');
hold on

quiver(X,Y,XMinusx(1:stride:end,1:stride:end),...
           XMinusy(1:stride:end,1:stride:end),2)
title('c^-');

subplot(2,2,3)
set(gca,'FontSize',16)
imagesc([0 Lx],[0 Ly], LambdaPlus)
colorbar
caxis([0 max(max(LambdaPlus))]);
set(gca,'YDir','normal');
title('\lambda^+');

subplot(2,2,4)
set(gca,'FontSize',16)
imagesc([0 Lx],[0 Ly], LambdaMinus)
colorbar
caxis([0 max(max(LambdaMinus))]);
set(gca,'YDir','normal');
title('\lambda^-');

