%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 02-May-2014 12:50:47
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [colors]=rainbow(R, G, B, l, L)
% R,G,B -> red,green,blue intensity (0:1)
% l -> loop variable
% L -> size of loop
function [colors]=rainbow(R, G, B, l, L)
om=2*pi/L;
red=R*sin(l*om-2*pi*(0/3));
green=G*sin(l*om-2*pi*(1/3));
blue=B*sin(l*om-2*pi*(2/3));
colors=[red green blue]/2+.5;
end
