%function [colors]=rainbow(R, G, B, l, L)
% om=2*pi/L;
% 
% 
% red=R*sin(l*om-2*pi*(0/3));
% green=G*sin(l*om-2*pi*(1/3));
% blue=B*sin(l*om-2*pi*(2/3));
% 
% colors=[red green blue]/2+.5;
function [colors]=rainbow(R, G, B, l, L)



om=2*pi/L;


red=R*sin(l*om-2*pi*(0/3));
green=G*sin(l*om-2*pi*(1/3));
blue=B*sin(l*om-2*pi*(2/3));

colors=[red green blue]/2+.5;

