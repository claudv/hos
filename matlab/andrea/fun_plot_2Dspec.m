function funout=fun_plot_2Dspec(date,datevec,specall,freq,dir)

% Plots the directional wave spectra for the timestep itime
% Returns the significant wave height for itime
%--------------------------------------------------------------------------
% INPUTS:
%  date: date string for the plot ('yyyymmdd HH:MM:SS')
%--- outputs of fun_2Dspec_read:
%  datevec: vec with the time stamps in num format (datenum)
%  specall: vec spectral components (dir,freq), for all time steps
%  freq: vec frequencies (in Hz)
%  dir:  vec directions (in radians)
%--------------------------------------------------------------------------
% the directions are given in the oceanographic convention 
%:(0 going to N, 90 going to E)
% stored in the cartesian 
%:(from 90->0->270->180->90)
% Note 1: this convention is used for makeboundary.f90 (conversion from
% .spec ASCII files to nest.ww3 binary files)
% Note 2: this should work for files with the directions in any order, as
% long as the ocenographic convention is respected -- NEEDS TO BE TESTED!!!
%--------------------------------------------------------------------------
% OUTPUTS
% significant wave height
%--------------------------------------------------------------------------
% Revision history
% May 22: sorting for both polar and cartesian, wrapping for polar
% July 24: reads spectrum for each time step; converts data
% July 25: converted to function
% September 19: reads in fun_2Dspec_read outputs rather than files
%--------------------------------------------------------------------------

cont1=0.1;maxs=100;
step=(maxs-cont1)/200;


%compute the number of dir, freq, timesteps
ndir=length(dir)
nfreq=length(freq)

ntime=length(datevec)

numDate=datenum(date)
itime=find(datevec==numDate);

%extract part of the spectrum for itime
spec=specall((itime-1)*nfreq+1:itime*nfreq,:);


ptitle=datestr([datevec(itime)]);
% sort directions in ascending order and then the spectral components
[dir_sort,dir_idx]=sort(dir');
spec_sort=zeros(nfreq,ndir);
for i=1:ndir
    spec_sort(:,i)=spec(:,dir_idx(i));
end


%-------------------------------------------------------------------
% Create polar data
dir_wrap=[dir_sort dir_sort(1)];
spec_wrap=[spec_sort spec_sort(:,1)];
[r,t] = meshgrid(freq,dir_wrap);
z = spec_wrap';
figure
set(gcf,'renderer','zbuffer')
contourf(freq',dir_wrap*180/pi,spec_wrap',[cont1:step:maxs],'EdgeColor','none');
xlim([freq(1) freq(end)]);
ylim([0 360])
caxis([cont1 maxs])
colorbar
%shading interp
title([ptitle,' ',': 2D wave spectra (m^2/Hz/rad)']);
xlabel('Freq (Hz)')
ylabel('Dir (deg)')
size(z);
% Convert to Cartesian
x = r.*sin(t);
y = r.*cos(t);

figure
h = polar(x,y);
hold on;
contourf(x,y,z,[cont1:step:maxs],'EdgeColor','None');

% Hide the POLAR function data and leave annotations
set(h,'Visible','off')
% Turn off axes and set square aspect ratio
axis off
axis image

%Removing the label
set(findall(gcf, 'String', '30', '-or','String','60','-or','String','90') ,'String', ' ');
set(findall(gcf, 'String', '330', '-or','String','360','-or','String','270') ,'String', ' ');
set(findall(gcf, 'String', '120', '-or','String','150','-or','String','180') ,'String', ' ');
set(findall(gcf, 'String', '210', '-or','String','240','-or','String','300') ,'String', ' ');

%Altering the label
set(findall(gcf, 'String', '0'),'String', ' 90'); 
set(findall(gcf, 'String', '90'),'String', ' 0');
colorbar;
caxis([cont1 maxs])

title([ptitle,' ',': 2D wave spectra (m^2/Hz/rad)']);


%--------------------------------------------------------------------
%--------------------------------------------------------------------
% COMPUTE INTEGRAL QUANTITIES


% periodic wrap

spec_periodic=[spec_sort spec_sort(:,1)];


hu=log(1.1);
hdir=2*pi/ndir;

integrant=spec_periodic;

for i=1:nfreq
integrant(i,:)= spec_periodic(i,:)*freq(i);
end

m0=0;

% finite part of the integral

for i=2:nfreq-1
    m0=m0+2.0*integrant(i,1)+2.0*integrant(i,ndir+1);
    for j=2:ndir
        m0=m0+4.0*integrant(i,j);
    end
end

for j=2:ndir
    m0=m0+2.0*integrant(1,j)+2.0*integrant(nfreq,j);
end

m0=m0+integrant(1,1)+integrant(1,ndir+1);
m0=m0+integrant(nfreq,1)+integrant(nfreq,ndir+1);

m0=m0*hdir*hu/4.0;

% spectral tail

tail=0.0;

for i=1:ndir+1
    tail=tail+spec_periodic(nfreq,i);
end

tail=hdir*(tail-0.5*spec_periodic(nfreq,1)-0.5*spec_periodic(nfreq,ndir+1));
tail=tail*0.25*freq(nfreq);

HS1=sqrt(m0+tail)*4;

funout.HS=HS1;

