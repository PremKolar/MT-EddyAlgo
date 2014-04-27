%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig(rez,xdim,ydim,tit)
	%% set up figure
	resolution=get(0,'ScreenPixelsPerInch');
	xdim=xdim*rez/resolution;
	ydim=ydim*rez/resolution;
	set(gcf,'paperunits','inch','papersize',[xdim ydim]/rez,'paperposition',[0 0 [xdim ydim]/rez]);
	mkdirp('~/FIGS/')
	fname=['~/FIGS/',tit,'.fig'];
	fnamepng=['~/FIGS/',tit,'.png'];
	saveas(gcf,fname);
	try
	saveas(gcf,fnamepng);
	end
	close(gcf);
end

