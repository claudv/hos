function [ CPlus, CMinus, LambdaPlus, LambdaMinus, XPlusx, XPlusy, XMinusx, XMinusy ] = LocalPhaseSpeed( eta, eta_t, W, Kx, Ky )

    [KX, KY] = meshgrid(Kx,Ky);
    heta = fftshift(fft2(eta));
    
    % eta_x --------------------%
    heta_x = 1i*KX.*heta;
    heta_x = ifftshift(heta_x);
    eta_x = real(ifft2(heta_x));

    % eta_y --------------------%
    heta_y = 1i*KY.*heta;
    heta_y = ifftshift(heta_y);
    eta_y = real(ifft2(heta_y));
    
    % W^2 ----------------------%
    WW = W.*W;
    hWW = fft2(WW);

    % fft2(eta_x.*eta_y) -------------------%
    heta_x_times_eta_y = fft2(eta_x.*eta_y);
    
    % fft2(eta_x.*eta_x) -------------------%
    heta_x_times_eta_x = fft2(eta_x.^2);
    
    % fft2(eta_y.*eta_y) -------------------%
    heta_y_times_eta_y = fft2(eta_y.^2);
    
    % fft2(eta_x.*eta_t) -------------------%
    heta_x_times_eta_t = fft2(eta_x.*eta_t);
    
    % fft2(eta_y.*eta_t) -------------------%
    heta_y_times_eta_t = fft2(eta_y.*eta_t);
    
    
    
    Pxy = real(ifft2(hWW.*heta_x_times_eta_y));
    Pxx = real(ifft2(hWW.*heta_x_times_eta_x));
    Pyy = real(ifft2(hWW.*heta_y_times_eta_y));
    Pxt = real(ifft2(hWW.*heta_x_times_eta_t));
    Pyt = real(ifft2(hWW.*heta_y_times_eta_t));
    
    
    % Eigenvalues ---------------------------------------------%
    LambdaPlus  = 0.5*(Pxx + Pyy) + sqrt( 0.25*(Pxx + Pyy).^2 - Pxx.*Pyy + Pxy.^2 );
    LambdaMinus = 0.5*(Pxx + Pyy) - sqrt( 0.25*(Pxx + Pyy).^2 - Pxx.*Pyy + Pxy.^2 );
    
    % Eigenvectors (principal directions) ---------------------%
    XPlusx = -Pxy./(Pxx-LambdaPlus);
    XPlusy = ones(size(XPlusx));
    Norm = abs(XPlusx + 1i*XPlusy);
    XPlusx = XPlusx./Norm;
    XPlusy = XPlusy./Norm;
    
    XMinusy = -Pxy./(Pyy-LambdaMinus);
    XMinusx = ones(size(XMinusy));
    Norm = abs(XMinusx + 1i*XMinusy);
    XMinusx = XMinusx./Norm;
    XMinusy = XMinusy./Norm;
    
    % Phase speed components along the principal directions ---%
    CPlus = -(XPlusx.*Pxt + XPlusy.*Pyt)./LambdaPlus; 
    CMinus = -(XMinusx.*Pxt + XMinusy.*Pyt)./LambdaMinus;
    
    
    Sign = sign(CPlus);
    XPlusx = Sign.*XPlusx;
    XPlusy = Sign.*XPlusy;
    CPlus = Sign.*CPlus;
    
    Sign = sign(CMinus);
    XMinusx = Sign.*XMinusx;
    XMinusy = Sign.*XMinusy;
    CMinus = Sign.*CMinus;
    
    % rescale eigenvalues --------------------------------------%
    %LambdaPlus = LambdaPlus/numel(Pxx);
    %LambdaMinus = LambdaMinus/numel(Pxx);
    
end

