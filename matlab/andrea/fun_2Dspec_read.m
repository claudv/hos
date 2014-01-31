function funout=fun_2Dspec_read(specname,nPoints,indexPoint)
%----------------------
% INPUTS:
% Reads a spec file (wavewatch form: specname
% the file can contain multiple points nPoints
% returns only spec data for point indexPoints
%------------------------
% OUTPUTS:
% spec(nfreq*ntime,ndir) spectral comp for all time steps
% date(ntime) a vector with time in numerical format
% freq(nfreq) freq
% dir(ndir) dir

fid=fopen(specname,'r');

% get frequencies and directions
out =fscanf(fid,'%*s %*s %*s %d %d',[1,6]);
frewind(fid);
fgetl(fid);
 
 nfreq=out(1);
 ndir=out(2);
 
 disp(sprintf('nfreq %d',nfreq))
 disp(sprintf('ndir  %d',ndir))
 
[freq]=fscanf(fid,'%e',[nfreq,1]);
[dir]=fscanf(fid,'%e',[ndir,1]);

%get the first date
fgetl(fid); %skip a line
tline=fgetl(fid);
spec=[];
date=[];

% scrol through the file for each date
while(ischar(tline))
   % disp('Getting data on ...');
   % disp(tline);
    dateN=datenum(tline,'yyyymmdd HHMMSS');
    date=[date;dateN];
    for i=1:nPoints % scroll through points
        tline=fgetl(fid);
        part=fscanf(fid,'%e',[nfreq,ndir]);
        if i==indexPoint
          % disp('For...');
          % disp(tline);
         spec=[spec; part];
        end
        fgetl(fid); %skip
    end
    tline=fgetl(fid);
      
end
fclose(fid);
funout.SPEC=spec;
funout.DATE=date;
funout.FREQ=freq;
funout.DIR=dir;
