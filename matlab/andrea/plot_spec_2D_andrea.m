%close all;
clear all;
clc;

date='9 November 2007 00:00:00';
numdate=datenum(date);


specdata=fun_2Dspec_read('andrea.spec',1,1);
specall=specdata.SPEC;
freq=specdata.FREQ;
dir=specdata.DIR;
nfreq=length(freq);
ndir=length(dir);
datevec=specdata.DATE;

fun_plot_2Dspec(date,datevec,specall,freq,dir)
%title('Var density spectra m^2/Hz/rad')


%% plot hs
hsfile='mean_wave_parameters_ANDREA.nc'
lat=ncread(hsfile,'g0_lat_1');
lon=ncread(hsfile,'g0_lon_2');
hs=ncread(hsfile,'SWH_GDS0_MSL');
time=ncread(hsfile,'initial_time0_hours');
t0=datenum('1 January 1800');
time=time/24+t0;

index_date=find(time==numdate);

figure
set(gcf,'renderer','zbuffer')
pcolor(lon',lat,hs(:,:,index_date)')
shading interp
colorbar
title(['hs (m) ' datestr(numdate)])
hold on
plot(3.2,56.5,'ro','MarkerFaceColor','black')
