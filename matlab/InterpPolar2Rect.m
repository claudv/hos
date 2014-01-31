%   Data interpolation from polar to rectangular grid.
%
%   [ xr, yr, vr] = InterpPolar2Rect( x, y, v );
%
%   Inputs ( length(x) = length(y) = length(v) )
%   x: vector of horizontal coordinates of the nodels of the polar grid.
%   y: vector of vertical coordinates of the nodes of the polar grid.
%   v: data on the polar grid.
%
%   Outputs
%   xr: vector of the x coordinates of the rectangular grid.
%   yr: vector of the y coordinates of the rectangular grid.
%   vr: data interpolated on the rectangular grid.

function [ xr, yr, vr] = InterpPolar2Rect( x, y, v )


    vq = griddata(x,y,v,xq,yq,'cubic');

    
    
    
end

