%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig(outdir,rez,xdim,ydim,tit)
    disp('yo')
	fname=[outdir,tit];
     mkdirp([outdir,'old']);	
% 	junkdir=[outdir,'old/', 'movedOn' datestr(now,'mmdd') '/'];
    mkdirp(outdir);	
%       mkdirp(junkdir);
%       system(['mv ' outdir '*pdf ' outdir '*png ' junkdir])
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
% 		disp(['printing' fname ' to ps']);
% 		print(gcf, '-dpsc2', [fname,'.ps'])
		eval(['print ',[fname,'.eps'] , ' -f -r',num2str(rez),' -depsc ;'])
		system(['epstopdf ' fname '.eps']);
% 		system(['pdfcrop ' fname '.pdf ' fname 'Crpd.pdf']);
		system(['rm ' fname '.eps']);
	end
	close(gcf);
end

