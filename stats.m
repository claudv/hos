%script for evaluating various statistics of HOS output data

close all

SimFolder = '~/rw/hos_data/1';

N = 23;



kurt  = zeros(1,N);
skew  = zeros(1,N);
t     = zeros(1,N);
edge  = -5:0.01:5;
gauss = exp(-0.5*(edge.^2));

for nfield=[0:1:N-1];
    
    time = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/time');
    eta  = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/eta');
    phi  = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/phi');
    Lx   = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Lx');
    Ly   = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ly');
    Nx   = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Nx');
    Ny   = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/Ny');
    g    = h5read([SimFolder,'/data',num2str(nfield),'.1.h5'],'/g');
   
    Kx = [-Nx/2:1:Nx/2-1]'*2*pi/Lx;
    Ky = [-Ny/2:1:Ny/2-1]'*2*pi/Ly;
    x  = Lx*(0:Nx-1)/Nx;
    y  = Ly*(0:Ny-1)/Ny;
    
    eta_vec = reshape(eta,1,Nx*Ny);
    
    i    = nfield+1;  
    t(i) = i; 
    
    mu3     = mean((eta_vec - mean(eta_vec)).^3);
    mu4     = mean((eta_vec - mean(eta_vec)).^4);
    skew(i) = mu3/std(eta_vec)^3;
    kurt(i) = mu4/std(eta_vec)^4;
    
        %normalisation
    eta_vec = (eta_vec - mean(eta_vec))/std(eta_vec);
    %edge    = min(eta_vec):0.01:max(eta_vec);
    pdf     = histc(eta_vec,edge);
    norm_g  = max(pdf)*gauss;
    
    figure(1)
   
    
    semilogy(edge,norm_g,'r','LineWidth',2)
    hold on
    semilogy(edge,pdf,'b')
    axis([-5 5 0.5 10^4])
    annotation('textbox',[0.63 0.8 0.25 0.10],'String',...
        {['skewness = ' num2str(skew(i))],['kurtosis = ' num2str(kurt(i))]},...
        'Fontsize', 12, 'Fontname', 'Arial', 'LineWidth', 2,...
        'BackgroundColor', [0.9 0.9 0.9]);
    
    title ('$\eta$ probabilty density function','interpreter','latex','fontsize',22);
    legend('Normal Distribution','\eta PDF','location','NorthWest')
    
    hold off
    pause();
    
end

figure(2)
plot(t*60,kurt,'bd-');
title ('Kurtosis evolution','interpreter','latex','fontsize',22);
xlabel ('time (s)','interpreter','latex','fontsize',18);
ylabel ('kurtosis','interpreter','latex','fontsize',18);

figure(3)
plot(t*60,skew,'ro-');
title ('Skewness evolution','interpreter','latex','fontsize',22);
xlabel ('time (s)','interpreter','latex','fontsize',18);
ylabel ('skewness','interpreter','latex','fontsize',18);

