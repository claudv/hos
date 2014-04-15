heta_r=h5read('data_extra0.1.h5','/heta_r');
heta_i=h5read('data_extra0.1.h5','/heta_i');
heta=heta_r + i*heta_i;

hetalow = flipud(heta(2:end,:));
hetalow(:,2:end)=fliplr(hetalow(:,2:end));
hetalow=conj(hetalow);


hetaall = [heta; zeros(1,size(heta,2)); hetalow];

hetaall =  hetaall*size(hetaall,1)*size(hetaall,2);

eta_small = real(ifft2(hetaall));

eta=h5read('data0.1.h5','/eta');