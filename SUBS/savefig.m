%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig(outdir,rez,xdim,ydim,tit,debug,vectorGr)
    if nargin<6,debug=false;end
    if debug
    disp('yo')
    return
	 end
	 if nargin < 7
		 vectorGr='dpng';
	 end
	fname=[outdir,tit];
    mkdirp(outdir);	
	if rez==42
		%% quick hack
		disp(['quick print to ' [fname,'.png']])
		print(gcf, '-dpng', [fname,'.png'])
			
	else
		%% set up figure
		set(gcf,'renderer','opengl');
		resolution=get(0,'ScreenPixelsPerInch');
		xdim=xdim*rez/resolution;
		ydim=ydim*rez/resolution;
		set(gcf,'paperunits','inch','papersize',[xdim ydim]/rez,'paperposition',[0 0 [xdim ydim]/rez]);
		if strcmp(vectorGr,'dpdf')
		eval(['print ',[fname,'.eps'] , ' -f -r',num2str(rez),' -depsc ;'])
		system(['epstopdf ' fname '.eps']);
		system(['rm ' fname '.eps']);
		else
			eval(['print ',[fname,vectorGr(2:end)] , ' -f -r',num2str(rez),' -',vectorGr,';'])
		end
	end
	close(gcf);
end

