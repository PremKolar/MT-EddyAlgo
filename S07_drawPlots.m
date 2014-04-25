%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Apr-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S07_drawPlots
	%% init
	close all
	DD=initialise('cuts');
	%%
	%DD.threads.lims_eddies=thread_distro(DD.threads.num,numel(DD.path.tracks.files));
	%%
	init_threads(DD.threads.num);
	maps=load([DD.path.analyzed.name, 'maps.mat']);
	la=maps.Cycs.GLA;
	lo=maps.Cycs.GLO;
	tracks=load([DD.path.analyzed.name, 'tracks.mat']);
	%%
	
	ticks.y=[-80 -60 -40 -20 0 20 40 60 80]';
	ticks.x=linspace(-180,180,9);
	ticks.age=[1,DD.time.span,10];
	ticks.isoper=[DD.thresh.shape.iq,1,10];
	ticks.radius=[1,200,10];
	ticks.amp=[0,30,10];
	ticks.visits=[0,max([maps.AntiCycs.visitsSingleEddy(:); maps.Cycs.visitsSingleEddy(:)]),10];
	ticks.dist=[min([maps.AntiCycs.dist.zonal.fromBirth.mean(:);...
		maps.Cycs.dist.zonal.fromBirth.mean(:)])/1000,...
		max([maps.AntiCycs.dist.zonal.fromBirth.mean(:); maps.Cycs.dist.zonal.fromBirth.mean(:)])/1000,10];
	ticks.dist=[-80;60;8];
	ticks.disttot=[0;200;10];
	ticks.vel=[-50;20;15];
	
	ticks.axis=[DD.map.geo.west DD.map.geo.east DD.map.geo.south DD.map.geo.north	];
	
	
	
	
	%spmd
		%%
	%	if labindex==2
			trackPlots(DD,ticks,tracks)
	%	end
	%	if labindex==1
			mapstuff(maps,DD,ticks,lo,la)
	%	end
	%end
		%% update infofile
		%	save_info(DD)

end

function mapstuff(maps,DD,ticks,lo,la)
	for sense=fieldnames(maps)';sen=sense{1};
		maps.(sen).age.logmean=log(maps.(sen).age.mean);
		%%
		figure
		VV=maps.(sen).age.logmean;
		pcolor(lo,la,VV);shading flat
		%caxis([0 max(VV(:))]);
		drawcoast;
		decorate('age',ticks,DD,sen,'age','d',1,1);
		%%
		figure
		VV=maps.(sen).visitsSingleEddy;
		pcolor(lo,la,VV);shading flat
		%caxis([0 nanmax(VV(:))]);
		drawcoast;
		decorate('visits',ticks,DD,sen,'Visits','',0,1);
		%%
		figure
		VV=maps.(sen).dist.zonal.fromBirth.mean/1000;
		pcolor(lo,la,VV);shading flat
		%caxis([0 nanmax(VV(:))]);
		drawcoast;
		cb=decorate('dist',ticks,DD,sen,'Distance from Birth','km',0,1);
		doublemap(cb,winter,hot,[1 1 1])
		%%
		figure
		VV=maps.(sen).dist.zonal.tillDeath.mean/1000;
		pcolor(lo,la,VV);shading flat
		%caxis([0 nanmax(VV(:))]);
		drawcoast;
		cb=decorate('dist',ticks,DD,sen,'Distance till Death','km',0,1);
		doublemap(cb,winter,hot,[1 1 1])
		%%
		figure
		VV=maps.(sen).dist.traj.fromBirth.mean/1000;
		pcolor(lo,la,VV);shading flat
		%caxis([0 nanmax(VV(:))]);
		drawcoast;
		decorate('disttot',ticks,DD,sen,'Total distance travelled since birth','km',0,1);
		%%
		figure
		VV=maps.(sen).dist.traj.tillDeath.mean/1000;
		pcolor(lo,la,VV);shading flat
		%caxis([0 nanmax(VV(:))]);
		drawcoast;
		decorate('disttot',ticks,DD,sen,'Total distance to be travelled till death','km',0,1);
		%%
		figure
		VV=maps.(sen).vel.zonal.mean*100;
		pcolor(lo,la,VV);shading flat
		%caxis([0 nanmax(VV(:))]);
		drawcoast;
		cb=decorate('vel',ticks,DD,sen,'Zonal velocity','cm/s',0,1);
		doublemap(cb,jet,bone,[1 .95 .95])
		%%
		figure
		VV=maps.(sen).radius.mean.mean/1000;
		pcolor(lo,la,VV);shading flat
		%caxis([0 nanmax(VV(:))]);
		drawcoast;
		decorate('radius',ticks,DD,sen,'Radius','km',0,1);
		
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
function trackPlots(DD,ticks,tracks)
	%
% 	figure
% 	drawColorLine(tracks.Cycs,'age',240,1) ;
% 	decorate('age',ticks,DD,'Cyclones','age','d',1,1)
	%%
	figure
	[ACyc.Fig,ACyc.maxTick,ACyc.cmap]=drawColorLine(tracks.AntiCycs,'age',240,1) ;
	ACyc.mintick=1;
	decorate('age',ticks,DD,'Anti-Cyclones','age','d',1,1)
	%%
% 	figure
% 	[Cyc.Fig,Cyc.maxTick,Cyc.cmap]=drawColorLine(tracks.Cycs,'isoper',1,0) ;
% 	Cyc.mintick=DD.thresh.shape.iq;
% 	decorate('isoper',ticks,DD,'Cyclones','IQ',' ',0,100)
	%%
	figure
	[ACyc.Fig,ACyc.maxTick,ACyc.cmap]=drawColorLine(tracks.AntiCycs,'isoper',1,0) ;
	ACyc.mintick=DD.thresh.shape.iq;
	decorate('isoper',ticks,DD,'Anti-Cyclones','IQ',' ',0,100)
	%%
% 	figure
% 	[Cyc.Fig,Cyc.maxTick,Cyc.cmap]=drawColorLine(tracks.Cycs,'radiusmean',max(cat(2,tracks.Cycs.radiusmean)),0) ;
% 	Cyc.mintick=1000;
% 	decorate('radius',ticks,DD,'Cyclones','Radius','km',0,1)
	%%
	figure
	[ACyc.Fig,ACyc.maxTick,ACyc.cmap]=drawColorLine(tracks.AntiCycs,'radiusmean',max(cat(2,tracks.AntiCycs.radiusmean)),0) ;
	ACyc.mintick=1000;
	decorate('radius',ticks,DD,'Anti-Cyclones','Radius', 'km',0,1)
	%%
% 	figure
% 	[Cyc.Fig,Cyc.maxTick,Cyc.cmap]=drawColorLine(tracks.Cycs,'peakampto_ellipse',max(cat(2,tracks.Cycs.peakampto_ellipse)),0) ;
% 	Cyc.mintick=0;
% 	decorate('amp',ticks,DD,'Cyclones','Amp to ellipse','cm',0,1)
	%%
	figure
	[ACyc.Fig,ACyc.maxTick,ACyc.cmap]=drawColorLine(tracks.AntiCycs,'peakampto_ellipse',max(cat(2,tracks.AntiCycs.peakampto_ellipse)),0) ;
	ACyc.mintick=0;
	decorate('amp',ticks,DD,'Anti-Cyclones','Amp to ellipse','cm',0,1)
	%%
% 	figure
% 	[Cyc.Fig,Cyc.maxTick,Cyc.cmap]=drawColorLine(tracks.Cycs,'peakampto_contour',max(cat(2,tracks.Cycs.peakampto_ellipse)),0) ;
% 	Cyc.mintick=0;
% 	decorate('amp',ticks,DD,'Cyclones','Amp to contour','cm',0,1)
	%%
	figure
	[ACyc.Fig,ACyc.maxTick,ACyc.cmap]=drawColorLine(tracks.AntiCycs,'peakampto_contour',max(cat(2,tracks.AntiCycs.peakampto_ellipse)),0) ;
	ACyc.mintick=0;
	decorate('amp',ticks,DD,'Anti-Cyclones','Amp to contour','cm',0,1)
	
	
	
	
	
end
function cb=decorate(field,ticks,DD,tit,tit2,unit,logornot,decim)
	set(gcf,'renderer','opengl')
	axis(ticks.axis);
	set(gca,'ytick',ticks.y);
	set(gca,'xtick',ticks.x) ;
	cb=colorbar;
	load coast;
	hold on; plot(long,lat,'LineWidth',0.5);
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
	
	printAtRes(400,500,450)
end
function [h,maxV,cmap]=drawColorLine(V,fieldName,maxV,logornot)
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
				%h(ii)=line([lo(ii) lo(ii+1)],[la(ii) la(ii+1)],'color',cm(:,ii),'LineWidth',0.5);
				line([lo(ii) lo(ii+1)],[la(ii) la(ii+1)],'color',cm(:,ii),'LineWidth',0.5);
			end
		end
	end
	h=nan;
	caxis([minV maxV])
end
