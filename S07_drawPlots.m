%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Apr-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all figs are printed to ~/PRINTS/  !
function S07_drawPlots
	%% init
	[DD,maps,tracks,lo,la]=inits;
	%%	set ticks here!
	ticks.rez=400;
	ticks.width=800;
	ticks.height=800;
	ticks.y=linspace(0,40,5);
	ticks.x=linspace(0,40,5);
	ticks.age=[1,DD.time.span/2,10];
	ticks.isoper=[DD.thresh.shape.iq,1,10];
	ticks.radius=[0,300,7];
	ticks.amp=[0,20,7];
	ticks.visits=[0,max([maps.AntiCycs.visitsSingleEddy(:); maps.Cycs.visitsSingleEddy(:)]),5];
	%ticks.dist=[min([maps.AntiCycs.dist.zonal.fromBirth.mean(:);...
	%		maps.Cycs.dist.zonal.fromBirth.mean(:)])/1000,...
	%		max([maps.AntiCycs.dist.zonal.fromBirth.mean(:); maps.Cycs.dist.zonal.fromBirth.mean(:)])/1000,10];
	ticks.dist=[-1200;300;8];
	ticks.disttot=[0;2600;8];
	ticks.vel=[-30;10;5];
	%ticks.axis=[DD.map.geo.west DD.map.geo.east DD.map.geo.south DD.map.geo.north	];
	ticks.axis=[0 46 0 37];
	%%
	main(DD,tracks,maps,lo,la,ticks);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,maps,tracks,lo,la]=inits
	DD=initialise('cuts');
	init_threads(2);
	maps=load([DD.path.analyzed.name, 'maps.mat']);
	la=maps.Cycs.GLA;
	lo=maps.Cycs.GLO;
	tracks=load([DD.path.analyzed.name, 'tracks.mat']);
end
function main(DD,tracks,maps,lo,la,ticks)
	spmd
		if labindex==2
			trackPlots(DD,ticks,tracks)
		end
		if labindex==1
			mapstuff(maps,DD,ticks,lo,la)
		end
	end
end
function mapstuff(maps,DD,ticks,lo,la)
	for sense=fieldnames(maps)';sen=sense{1};
		maps.(sen).age.logmean=log(maps.(sen).age.mean);
		%%
		%figure
		VV=maps.(sen).age.logmean;
		pcolor(lo,la,VV);shading flat
		decorate('age',ticks,DD,sen,'age','d',1,1);
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapAge'])
		%%
		%figure
		VV=maps.(sen).visitsSingleEddy;
		pcolor(lo,la,VV);shading flat
		decorate('visits',ticks,DD,sen,'Visits of unique eddy','',0,1);
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapVisits']);
		%%
		%figure
		VV=maps.(sen).dist.zonal.fromBirth.mean/1000;
		pcolor(lo,la,VV);shading flat
		cb=decorate('dist',ticks,DD,sen,'Distance from Birth','km',0,1);
		doublemap(cb,winter,autumn,[.9 1 .9])
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapDFB']);
		
		%%
		%figure
		VV=maps.(sen).dist.zonal.tillDeath.mean/1000;
		pcolor(lo,la,VV);shading flat
		cb=decorate('dist',ticks,DD,sen,'Distance till Death','km',0,1);
		doublemap(cb,winter,autumn,[.9 1 .9])
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapDTD']);
		%%
		%figure
		VV=maps.(sen).dist.traj.fromBirth.mean/1000;
		pcolor(lo,la,VV);shading flat
		decorate('disttot',ticks,DD,sen,'Total distance travelled since birth','km',0,1);
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapDTFB']);
		%%
		%figure
		VV=maps.(sen).dist.traj.tillDeath.mean/1000;
		pcolor(lo,la,VV);shading flat
		decorate('disttot',ticks,DD,sen,'Total distance to be travelled till death','km',0,1);
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapDTTD']);
		%%
		%figure
		VV=maps.(sen).vel.zonal.mean*100;
		pcolor(lo,la,VV);shading flat
		cb=decorate('vel',ticks,DD,sen,'Zonal velocity','cm/s',0,1);
		doublemap(cb,autumn,winter,[.9 1 .9])
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapVel']);
		%%
		%figure
		VV=maps.(sen).radius.mean.mean/1000;
		pcolor(lo,la,VV);shading flat
		decorate('radius',ticks,DD,sen,'Radius','km',0,1);
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapRad']);
		
	end
	
end
function trackPlots(DD,ticks,tracks)
	%%
	senses={'Cycs','AntiCycs'};
	%figure
	for sense=senses; sen=sense{1};
		field='age';
		drawColorLine(tracks.(sen),field,240,1) ;
		decorate(field,ticks,DD,sen,field,'d',1,1);
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_field']);
		%%
		%figure
		drawColorLine(tracks.(sen),'isoper',1,0) ;
		decorate('isoper',ticks,DD,sen,'IQ',' ',0,100)
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_IQ']);
		%%
		%figure
		drawColorLine(tracks.(sen),'radiusmean',max(cat(2,tracks.(sen).radiusmean)),0) ;
		decorate('radius',ticks,DD,sen,'Radius','km',0,1)
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_radius']);
		%%
		%figure
		drawColorLine(tracks.(sen),'peakampto_ellipse',max(cat(2,tracks.(sen).peakampto_ellipse)),0) ;
		decorate('amp',ticks,DD,sen,'Amp to ellipse','cm',0,1)
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_TrackPeakampto_ellipse']);
		%%
		%figure
		drawColorLine(tracks.(sen),'peakampto_contour',max(cat(2,tracks.(sen).peakampto_ellipse)),0) ;
		decorate('amp',ticks,DD,sen,'Amp to contour','cm',0,1)
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_TrackPeakampto_contour'])
	end
end

function doublemap(cb,cm1,cm2,centercol)
	%% get colorbardata
	zlim=(get(cb,'ylim'));
	ztick=(get(cb,'ytick'));
	zticklabel=(get(cb,'yticklabel'));
	%% resample to fit ticks
	cm1r=resample(cm1,abs(zlim(1)),1);
	cm2r=resample(cm2,zlim(2),1);
	CM=[cm1r;flipud(cm2r)];
	%% blend in the middle
	midfilt=linspace(-1,zlim(2)/abs(zlim(1)),length(CM));
	gp=repmat(gauspuls(midfilt,1,1,-1/5)',1,3);
	centercolvec=repmat(centercol,size(CM,1),1);
	CM=(1-gp).*CM + gp.*centercolvec;
	%% correct for round errors
	CM(CM<0)=0;
	CM(CM>1)=1;
	%% reset to old params
	colormap(CM);
	caxis(zlim);
	set(cb,'ytick',ztick);
	set(cb,'yticklabel',zticklabel);
end
function drawcoast
	load coast;
	hold on; plot(long,lat,'LineWidth',0.5);
end
function cb=decorate(field,ticks,DD,tit,tit2,unit,logornot,decim)
	set(gcf,'renderer','opengl')
	axis(ticks.axis);
	set(gca,'ytick',ticks.y);
	set(gca,'xtick',ticks.x) ;
	cb=colorbar;
	% 	load coast;
	hold on;
	% 	plot(long,lat,'LineWidth',0.5);
	if logornot
		zticks=linspace(log(ticks.(field)(1)),log(ticks.(field)(2)),ticks.(field)(3))';
		zticklabel=num2str(round(exp(zticks)));
	else
		zticks=linspace(ticks.(field)(1),ticks.(field)(2),ticks.(field)(3))';
		zticklabel=num2str(round(zticks*decim)/decim);
	end
	caxis([zticks(1) zticks(end)])
	set(cb,'ytick',zticks);
	set(cb,'yticklabel',zticklabel);
	title([tit,' - ',tit2,' [',unit,']'])
	xlabel(['Eddies that died younger ',num2str(DD.thresh.life),' days are excluded'])
end
function [maxV,cmap]=drawColorLine(V,fieldName,maxV,logornot)
	cmap=jet;% Generate range of color indices that map to cmap
	if logornot
		maxV=log(maxV);
		minV=1;
	else
		minV=min(cat(2,V.(fieldName)));
	end
	kk=linspace(minV,maxV,size(cmap,1));
	for ee=1:numel(V)
		VV=V(ee).(fieldName);
		if isempty(VV)
			continue
		end
		if logornot
			VV(VV==0)=1;
			cm = spline(kk,cmap',log(VV));  % Find interpolated colorvalue
		else
			cm = spline(kk,cmap',VV);       % Find interpolated colorvalue
		end
		cm(cm>1)=1;                        % Sometimes iterpolation gives values that are out of [0,1] range...
		cm(cm<0)=0;
		lo=V(ee).lon;
		la=V(ee).lat;
		for ii=1:length(la)-1
			if  abs(lo(ii+1)-lo(ii))<10 % avoid 0->360 jumps
				line([lo(ii) lo(ii+1)],[la(ii) la(ii+1)],'color',cm(:,ii),'LineWidth',0.5);
			end
		end
	end
	caxis([minV maxV])
end
