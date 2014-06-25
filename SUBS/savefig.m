%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig(outdir,rez,xdim,ydim,tit,debug,frmt,saveFig)
	if nargin<6,debug=false;end
	if debug,disp('yo');	return;	end
	if nargin < 7,	frmt='dpng';	end
	if nargin < 8,	saveFig=false;	end
	fname=[outdir,tit];
	mkdirp(outdir);
	if rez==42
		quickhack(fname)
	else
		%% set up figure
		setupfigure(rez,xdim,ydim)
		%% print
		printStuff(frmt,fname,rez,xdim,ydim)
	end
	if saveFig
		save(gcf,[fname,'.mat'])
	end
	close(gcf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printStuff(frmt,fname,rez,xdim,ydim)
	if strcmp(frmt,'dpdf')
		eval(['print ',[fname,'.eps'] , ' -f -r',num2str(rez),' -depsc ;'])
		system(['epstopdf ' fname '.eps']);
		system(['rm ' fname '.eps']);
	else
		fnfull=[fname,'.',frmt(2:end)];
		eval(['print ',fnfull , ' -f -r',num2str(rez),' -',frmt,';'])
% 		system(['convert -density ' num2str(rez) 'x' num2str(rez) ' -resize ' num2str(xdim) 'x' num2str(ydim) ' quality 100 ' fnfull ' ' fname '.pdf' ]);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setupfigure(rez,xdim,ydim)
	set(gcf,'renderer','opengl');
	resolution=get(0,'ScreenPixelsPerInch');
	xdim=xdim*rez/resolution;
	ydim=ydim*rez/resolution;
	set(gcf,'paperunits','inch','papersize',[xdim ydim]/rez,'paperposition',[0 0 [xdim ydim]/rez]);
	xa4=11.692*resolution;
	fsScaled=round(12/xa4*xdim)		;
	set(gca,'FontSize',fsScaled)
	set(findall(gcf,'type','text'),'FontSize',fsScaled)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quickhack(fname)
	disp(['quick print to ' [fname,'.png']])
	print(gcf, '-dpng', [fname,'.png'])
end



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	xoy=xdim/ydim;
%	ya4=xa4/xoy;
% set(gcf,'paperunits','inch','position',[0 0 [xdim ydim]]);
