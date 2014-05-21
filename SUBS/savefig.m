%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig(outdir,rez,xdim,ydim,tit,zbub)
%% set up figure
if nargin<6
    zbub=false;
end
set(gcf,'renderer','opengl');
resolution=get(0,'ScreenPixelsPerInch');
xdim=xdim*rez/resolution;
ydim=ydim*rez/resolution;
set(gcf,'paperunits','inch','papersize',[xdim ydim]/rez,'paperposition',[0 0 [xdim ydim]/rez]);
mkdirp(outdir);
<<<<<<< HEAD
fnamepng=[outdir,tit,'.png'];
fnameeps=[outdir,tit,'.eps'];
fnamepdf=[outdir,tit,'.pdf'];

if zbub    
    set(gcf,'renderer','zbuffer');
    eval(['print ', fnamepng, ' -f -r',num2str(rez),' -dpng ;'])
    system(['convert -quality 100 -density 30 ' fnamepng,' ',fnamepdf]);
else
    eval(['print ', fnameeps, ' -f -r',num2str(rez),' -depsc ;'])
    system(['epstopdf ' fnameeps]);
end
=======
fname=[outdir,tit];

% fnamepdf=[outdir,tit,'.pdf'];

% if zbub    
%     set(gcf,'renderer','zbuffer');
disp(['printing' fname ' to ps']);
    print(gcf, '-dpsc2', [fname,'.ps'])
     
eval(['print ',[fname,'.eps'] , ' -f -r',num2str(rez),' -depsc ;'])
 system(['epstopdf ' fname '.eps']);
 system(['pdfcrop ' fname '.pdf ' fname 'Crpd.pdf']);
 system(['rm ' fname '.eps']);
 system(['rm ' fname '.pdf']);
    
    
%     system(['convert -quality 100 -density 30 ' fnamepng,' ',fnamepdf]);
% else
%    
%    
% end
>>>>>>> master
close(gcf);
end

