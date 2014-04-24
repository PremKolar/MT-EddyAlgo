function printAtRes(rez,xdim,ydim)
	
	%% set up figure
	resolution=get(0,'ScreenPixelsPerInch');
	xdim=xdim*rez/resolution;
	ydim=ydim*rez/resolution;
	set(gcf,'paperunits','inch','papersize',[xdim ydim]/rez,'paperposition',[0 0 [xdim ydim]/rez]);
	
	%% print
	mkdirp('~/PRINTS/')
	fname=['~/PRINTS/print_',datestr(now,'mmddHHMMSS')];
	fnamepng=[fname,'.png'];
	fnamepdf=[fname,'.pdf'];
	print(gcf, '-dpng',['-r',num2str(rez)],fnamepng );
	print(gcf, '-dpdf',['-r',num2str(rez)],fnamepdf );
	
end

