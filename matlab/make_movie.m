clf
axes('Position',[0 0 1 1])
%movie(Movie)

writerObj = VideoWriter('movie.1','Motion JPEG AVI');
%writerObj.CompressionRatio = 30;	% Only with Motion JPEG 2000 profiles
writerObj.Quality = 75;			% Only with MPEG-4 or Motion JPEG AVI profile
writerObj.FrameRate = 10;		% Frames per second
open(writerObj);
MovieSmall=Movie(:,1:2:end);
writeVideo(writerObj,MovieSmall);
close(writerObj);