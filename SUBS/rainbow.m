function [colors]=rainbow(R, G, B, l, L)
    om=2*pi/L;
    red=sin(l*om-2*pi*(0/3))';
    green=sin(l*om-2*pi*(1/3))';
    blue=sin(l*om-2*pi*(2/3))';
    colors=([red green blue]/2+.5).*repmat([R G B],numel(l),1);
end
