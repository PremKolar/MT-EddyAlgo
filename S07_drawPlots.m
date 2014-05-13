%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Apr-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all figs are saved to ~/FIGS/  !
function S07_drawPlots
	%% init
	[DD,maps,tracks,vecs,lo,la]=inits;
	%%	set ticks here!
	ticks.rez=300;
	ticks.width=300;
	ticks.height=300;
	% 	ticks.y= linspace(20,60,2);
	% 	ticks.x= linspace(-90,-50,2);
	ticks.y= 0;
	ticks.x= 0;
	ticks.age=[1,DD.time.span/2,10];
	ticks.isoper=[DD.thresh.shape.iq,1,10];
	ticks.radius=[50,150,6];
	ticks.amp=[1,20,7];
	%ticks.visits=[0,max([maps.AntiCycs.visitsSingleEddy(:); maps.Cycs.visitsSingleEddy(:)]),5];
	ticks.visits=[1,20,11];
	ticks.dist=[-1400;300;8];
	ticks.disttot=[1;2000;14];
	ticks.vel=[-30;20;6];
	ticks.axis=[DD.map.geo.west DD.map.geo.east DD.map.geo.south DD.map.geo.north	];
	ticks.lat=[ticks.axis(3:4),5];
	
	
	
	%%
	main(DD,tracks,vecs,maps,lo,la,ticks);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,maps,tracks,vecs,lo,la]=inits
	% DD=initialise('cuts');
	DD=initialise();
	% init_threads(3);
	maps=load([DD.path.analyzed.name, 'maps.mat']);
	la=maps.Cycs.GLA;
	lo=maps.Cycs.GLO;
	
	%% collect tracks
	root=DD.path.analyzedTracks.AC.name;
	ACs={DD.path.analyzedTracks.AC.files.name};
	for ff=1:numel(ACs)
		tracks.AntiCycs(ff)=load([root ACs{ff}]);
	end
	root=DD.path.analyzedTracks.C.name;
	Cs={DD.path.analyzedTracks.C.files.name};
	for ff=1:numel(Cs)
		tracks.Cycs(ff)=load([root Cs{ff}]);
	end
	
	%% get vectors
	vecs=load([DD.path.analyzed.name, 'vecs.mat']);
end

function main(DD,tracks,vecs,maps,lo,la,ticks)
	spmd
		if labindex==1
			histstuff(vecs,DD,ticks)
		end
		if labindex==2
			trackPlots(DD,ticks,tracks)
		end
		if labindex==3
			mapstuff(maps,vecs,DD,ticks,lo,la)
		end
	end
end
function histstuff(vecs,DD,ticks)
for sense=fieldnames(vecs)';sen=sense{1};
    figure
    lat=vecs.(sen).lat;
    range.(sen).lat=round(min(vecs.(sen).lat)):round(max(vecs.(sen).lat));
    hc.(sen).lat=histc(lat,range.(sen).lat);
    semilogy(range.(sen).lat,hc.(sen).lat)
    tit=['number of ',sen,' per 1° lat']
    title([tit])
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_latNum']);
    %%
    figure
    age=vecs.(sen).age;
    range.(sen).age=round(min(vecs.(sen).age)):round(max(vecs.(sen).age));
    hc.(sen).age=histc(age,range.(sen).age);
    semilogy(range.(sen).age,hc.(sen).age)
    tit=['number of ',sen,' per age'];
    unit='d';
    title([tit,' [',unit,']'])
    hc.(sen).cum.age=fliplr(cumsum(fliplr(hc.(sen).age)));
    semilogy(range.(sen).age,hc.(sen).cum.age)
    tit=['upper tail cumulative of ',sen,' per age'];
    title([tit,' [',unit,']'])
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_ageUpCum']);
end
%%
figure
len=min([length(hc.Cycs.lat),length(hc.AntiCycs.lat)])
hc.rat.lat=hc.Cycs.lat(1:len)./hc.AntiCycs.lat(1:len);
hc.rat.lat(isinf(hc.rat.lat))=nan;
hc.rat.lat((hc.rat.lat==0))=nan;
totnum=hc.Cycs.lat(1:len) + hc.AntiCycs.lat(1:len);
x=(range.Cycs.lat(1:len));
y=log10(hc.rat.lat);
x(totnum==0)=[];
y(totnum==0)=[];
totnum(totnum==0)=[];
colors = flipud(bone(max(totnum)));
cols=colors(totnum,:);
for ii = 1:length(x)
  a=bar(x(ii), y(ii));
  hold on
  set(a,'facecolor', cols(ii,:));
end

axis([range.Cycs.lat(1) range.Cycs.lat(end) -1 1]);
ticks.y=linspace(-1,1,3)
ticks.lab.y=num2cell(10.^(ticks.y));
set(gca,'ytick',ticks.y)
set(gca,'yticklabel',ticks.lab.y)
tit=['ratio of cyclones to anticyclones as function of latitude'];
title(tit)
cm=flipud(bone)
colormap(cm)
cb=colorbar;
    zticks=linspace(0,max(totnum),5)';
    zticklabel=num2str(round((zticks)));   
caxis([zticks(1) zticks(end)])
set(cb,'ytick',zticks);
set(cb,'yticklabel',zticklabel);



savefig(ticks.rez,ticks.width,ticks.height,[sen,'_RatLat']);
%%
figure
len=min([length(hc.AntiCycs.age), length(hc.Cycs.age)]);
hc.rat.age=hc.Cycs.age(1:len)./hc.AntiCycs.age(1:len);
hc.rat.age(isinf(hc.rat.age))=nan;
hc.rat.age((hc.rat.age==0))=nan;
totnum=hc.Cycs.age(1:len) + hc.AntiCycs.age(1:len);
x=sqrt(range.Cycs.age(1:len));
y=log10(hc.rat.age);
x(totnum<2)=[];
y(totnum<2)=[];
totnum(totnum<2)=[];
tnlog=log(totnum)
tmin=min(tnlog);
tmax=max(tnlog);
cmap=jet;
kk=(linspace(tmin,tmax,size(cmap,1)));

for ii = 1:length(x)
  a=bar(x(ii), y(ii),.1);
  hold on
    cm = spline(kk,cmap',(tnlog(ii)));  % Find interpolated colorvalue
  cm(cm<0)=0;
  cm(cm>1)=1;
    set(a,'facecolor', cm');
end



axis([min(x)-1 max(x) -1 1]);
xtck=get(gca,'xtick')
set(gca,'xticklabel',xtck.^2)
ticks.y=linspace(-1,1,3)
ticks.lab.y=num2cell(10.^(ticks.y));
set(gca,'ytick',ticks.y)
set(gca,'yticklabel',ticks.lab.y)
tit=['ratio of cyclones to anticyclones as function of age'];
xlab=['age [d] - gray scale indicates total number of eddies available'];
title(tit)
xlabel(xlab)
cm=(autumn)
colormap(cm)
cb=colorbar;
    zticks=linspace(0,max(exp(totnum)),5)';
    zticklabel=num2str(round((zticks)));   
caxis([zticks(1) zticks(end)])
set(cb,'ytick',zticks);
set(cb,'yticklabel',zticklabel);
savefig(ticks.rez,ticks.width,ticks.height,[sen,'_LatNum']);
%%
figure
hc.rat.cum.age=hc.Cycs.cum.age(1:len)./hc.AntiCycs.cum.age(1:len);
hc.rat.cum.age(isinf(hc.rat.cum.age))=nan;
plot(range.Cycs.age(1:len),log10(hc.rat.cum.age));
axis([range.Cycs.age(1) range.Cycs.age(len) -1 1]);
ticks.y=linspace(-1,1,3);
ticks.lab.y=num2cell(10.^(ticks.y));
set(gca,'ytick',ticks.y)
set(gca,'yticklabel',ticks.lab.y)
tit=['ratio of cyclones to anticyclones as upper tail cum. of age'];
title(tit)
savefig(ticks.rez,ticks.width,ticks.height,[sen,'_RatAge']);
end


function mapstuff(maps,vecs,DD,ticks,lo,la)
	for sense=fieldnames(maps)';sen=sense{1};
		
		figure
		b.la=vecs.(sen).birth.lat;
		b.lo=vecs.(sen).birth.lon;
		d.la=vecs.(sen).death.lat;
		d.lo=vecs.(sen).death.lon;
		plot(b.lo,b.la,'.r',d.lo,d.la,'.g','markersize',5)
		hold on
		drawcoast
		axis(ticks.axis);
		set(gca,'ytick',ticks.y);
		set(gca,'xtick',ticks.x) ;
		legend('births','deaths')
		savefig(ticks.rez,ticks.width,ticks.height,[sen,'_deathsBirths']);
		
		
		
		
		%%
		%figure
		maps.(sen).age.logmean=log(maps.(sen).age.mean);
		VV=maps.(sen).age.logmean;
		pcolor(lo,la,VV);shading flat
		decorate('age',ticks,DD,sen,'age','d',1,1);
		savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapAge'])
		%%
		%figure
		VV=maps.(sen).visits.single;
		VV(VV==0)=nan;
		pcolor(lo,la,log(VV));shading flat
		decorate('visits',ticks,DD,sen,'Visits of unique eddy','',1,1);
		savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapVisits']);
		%%
		%figure
		VV=maps.(sen).visits.birth;
		VV(VV==0)=nan;
		pcolor(lo,la,log(VV));shading flat
		decorate('visits',ticks,DD,sen,'Births','',1,1);
		savefig(ticks.rez,ticks.width,ticks.height,[sen,'_births']);
		%%
		%figure
		VV=maps.(sen).visits.birth;
		VV(VV==0)=nan;
		pcolor(lo,la,log(VV));shading flat
		decorate('visits',ticks,DD,sen,'Births','',1,1);
		savefig(ticks.rez,ticks.width,ticks.height,[sen,'_deaths']);
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
		VV=log(maps.(sen).dist.traj.fromBirth.mean/1000);
		pcolor(lo,la,VV);shading flat
		decorate('disttot',ticks,DD,sen,'Total distance travelled since birth','km',0,1);
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_MapDTFB']);
		%%
		%figure
		VV=log(maps.(sen).dist.traj.tillDeath.mean/1000);
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
		
		drawColorLinem(tracks.(sen),'lat','isoper') ;
		title([sen '- deflections'])
		axis([-1200 300 -300 300])
		set(gca,'ytick',[-200 0 200])
		set(gca,'xtick',[-1000 0 200])
		colorbar
		xlabel('latitude repr. by color; IQ repr. by thickness')
		axis equal
		
		savefig(ticks.rez,ticks.width,ticks.height,[sen,'_defletcs']);
		
		
		
		%%
		field='age';
		drawColorLine(tracks.(sen),field,240,1,0) ;
		decorate(field,ticks,DD,sen,field,'d',1,1);
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_field']);
		%%
		%figure
		drawColorLine(tracks.(sen),'isoper',1,0,0) ;
		decorate('isoper',ticks,DD,sen,'IQ',' ',0,100)
		savefig(ticks.rez,ticks.width,ticks.height,['Mad',sen,'_IQ']);
		%%
		%figure
		drawColorLine(tracks.(sen),'radiusmean',max(cat(2,tracks.(sen).radiusmean)),0,0) ;
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
function cb=decorate(field,ticks,DD,tit,tit2,unit,logornot,decim,coast)
	if nargin<9
		coast=true;
	end
	axis(ticks.axis);
	set(gca,'ytick',ticks.y);
	set(gca,'xtick',ticks.x) ;
	cb=colorbar;
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
	if coast
		drawcoast;
	end
end

function [maxV,cmap]=drawColorLinem(V,fieldName,fieldName2)
	cmap=jet;% Generate range of color indices that map to cmap
	minV=min(cat(2,V.(fieldName)));
	maxV=max(cat(2,V.(fieldName)));
	minIQ=min(cat(2,V.(fieldName2)));
	maxIQ=max(cat(2,V.(fieldName2)));
	
	iqiq=linspace(minIQ,maxIQ,10);
	
	kk=linspace(minV,maxV,size(cmap,1));
	
	for ee=1:numel(V)
		VV=V(ee).(fieldName);
		VViq=V(ee).(fieldName2);
		if isempty(VV)
			continue
		end
		cm = spline(kk,cmap',VV);       % Find interpolated colorvalue
		cm(cm>1)=1;                        % Sometimes iterpolation gives values that are out of [0,1] range...
		cm(cm<0)=0;
		
		iq = spline(iqiq,linspace(0.01,1.5,10),VViq);       % Find interpolated thickness
		
		lo=V(ee).lon;
		la=V(ee).lat;
		%% deg2km
		yy=[0 cumsum(deg2km(diff(la)))];
		xx=[0 cumsum(deg2km(diff(lo)).*cosd((la(1:end-1)+la(2:end))/2))];
		for ii=1:length(xx)-1
			if  abs(xx(ii+1)-xx(ii))<1000 % avoid 0->360 jumps
				line([xx(ii) xx(ii+1)],[yy(ii) yy(ii+1)],'color',cm(:,ii),'LineWidth',iq(ii));
			end
		end
	end
	caxis([minV maxV])
end

function [maxV,cmap]=drawColorLine(V,fieldName,maxV,logornot,zeroshift)
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
		if zeroshift
			lo=lo-lo(1);
			la=la-la(1);
		end
		
		for ii=1:length(la)-1
			if  abs(lo(ii+1)-lo(ii))<10 % avoid 0->360 jumps
				line([lo(ii) lo(ii+1)],[la(ii) la(ii+1)],'color',cm(:,ii),'LineWidth',0.8);
			end
		end
	end
	caxis([minV maxV])
end
