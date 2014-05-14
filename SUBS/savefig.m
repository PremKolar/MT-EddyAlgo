%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig(rez,xdim,ydim,tit,zbub)
%% set up figure
if nargin<5
    zbub=false;
end

outdir='~/FIGS/W2/';


set(gcf,'renderer','opengl');
resolution=get(0,'ScreenPixelsPerInch');
xdim=xdim*rez/resolution;
ydim=ydim*rez/resolution;
set(gcf,'paperunits','inch','papersize',[xdim ydim]/rez,'paperposition',[0 0 [xdim ydim]/rez]);
mkdirp(outdir)
fnamepng=[outdir,tit,'.png'];
fnameeps=[outdir,tit,'.eps'];
fnamepdf=[outdir,tit,'.pdf'];

if zbub    
    set(gcf,'renderer','zbuffer');
    eval(['print ', fnamepng, ' -f -r',num2str(rez),' -dpng'])
    system(['convert -quality 100 -density 30 ' fnamepng,' ',fnamepdf])
else
    eval(['print ', fnameeps, ' -f -r',num2str(rez),' -depsc'])
    system(['epstopdf ' fnameeps])
end


close(gcf);
end

