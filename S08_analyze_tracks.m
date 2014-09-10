%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S08_analyze_tracks
	%% init 
   
	DD=initialise([],mfilename);
   
	%%
	DD.threads.tracks=thread_distro(DD.threads.num,numel(DD.path.tracks.files));
	%%
 	[map,vecs,minMax,scat]=main(DD);
	%%
	seq_body(minMax,map,DD,vecs,scat);
	%%
	conclude(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [map,vecs,minMax,scat]=main(DD)
	if DD.debugmode
		[map,vecs,minMax,scat]=spmd_body(DD);
	else
		spmd(DD.threads.num)
			[map,vecs,minMax,scat]=spmd_body(DD);
          
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main functions
function [MAP,V,MinMax,scatterStuff]=spmd_body(DD)
	%% get stuff	
	[MAP,V,JJ,MinMax]=initAll(DD);
	%%
	T=disp_progress('init','analyzing tracks');
	
    scatterStuff=struct;
 nn=0;
    for jj=JJ;
		T=disp_progress('calc',T,numel(JJ),100);
		%% get track
		[TT]=getTrack(DD,jj); if isempty(TT),continue;end
		%% chelt over IQ
        nn_old=nn;
     nn=nn+numel(TT.eddy.trck);
     newIdx=nn_old+1:nn;
       la=num2cell(extractfield(cat(1,TT.eddy.trck.geo),'lat'))';
%       area=num2cell(extractfield(cat(1,TT.eddy.trck.area),'total'))';      
%       [CoIQ(newIdx,1).IQ]    =deal(TT.eddy.trck.isoper);    
%       [CoIQ(newIdx,1).chelt] =deal(TT.eddy.trck.chelt);      
%       [CoIQ(newIdx,1).lat]   =deal(la{:}); 
%       [CoIQ(newIdx,1).area]  =deal(area{:}); 
      
      
       [scatterStuff(newIdx,1).IQ]    =deal(TT.eddy.trck.isoper);      
      [scatterStuff(newIdx,1).lat]   =deal(la{:}); 
%       [CoIQ(newIdx,1).areaTotal]  =extractdeepfield(TT,'eddy.trck.area.total')';
   tmp=num2cell(extractdeepfield(TT,'eddy.trck.area.RadiusOverRossbyL')');
[scatterStuff(newIdx,1).RoL]  =tmp{:};
     ptmp =extractdeepfield(TT,'eddy.trck.peak');      
     tmp=num2cell(extractdeepfield((cat(1,ptmp{:})),'amp.to_ellipse')');    
      [scatterStuff(newIdx,1).amp2ellip]  =tmp{:};
    
      
      
      
      
        %% resort tracks for output
		[MinMax]=resortTracks(DD,MinMax,TT);
		%% mapstuff prep
		switch TT.sense
			case -1
				[MAP.AntiCycs,V.AntiCycs]=MeanStdStuff( TT.eddy,MAP.AntiCycs,V.AntiCycs,DD);
			case 1
				[MAP.Cycs,V.Cycs]=MeanStdStuff( TT.eddy,MAP.Cycs,V.Cycs,DD);
		end
    end
    scatterStuff=gop(@vertcat, scatterStuff);
	%% get global extrms
	MinMax=globalExtr(MinMax);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MAP,V]=MeanStdStuff(eddy,MAP,V,DD)
	[MAP.strctr, eddy]=TRstructure(MAP,eddy,DD);
	if isempty(eddy.trck),return;end % out of bounds
	[NEW.age]=TRage(MAP,eddy);
	[NEW.dist,eddy]=TRdist(MAP,eddy);
	NEW.vel=TRvel(MAP,eddy);
	NEW.radius=TRradius(MAP,eddy);
	NEW.amp=TRamp(MAP,eddy);
	[NEW.visits.all,NEW.visits.single]=TRvisits(MAP);
	MAP=comboMS(MAP,NEW,DD);
	[V]=getVecs(eddy,V);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ACs,Cs]=netVels(DD,map)
	velmean=reshape(extractdeepfield(load(DD.path.meanU.file),...
		'means.small.zonal'),DD.map.out.Y,DD.map.out.X);
	ACs=	map.AntiCycs.vel.zonal.mean -velmean;
	Cs =	map.Cycs.vel.zonal.mean		 -velmean;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq_body(minMax,map,DD,vecs,scat)
	[map,vecs]=mergeThreadData(minMax,map,DD,vecs);  %#ok<NASGU>
	%% get rossby radius
	if DD.switchs.RossbyStuff
		map.Rossby=loadRossby(DD);        
		%% build radius/rossbyRadius ratio
		map.AntiCycs.radius.toRo=map.AntiCycs.radius.mean.mean./map.Rossby.small.radius;
        map.AntiCycs.radius.to2Ro=map.AntiCycs.radius.mean.mean./(2*map.Rossby.small.radius);
		map.Cycs.radius.toRo=map.Cycs.radius.mean.mean./map.Rossby.small.radius;
	end
	%% build zonal means
	map.zonMean=zonmeans(map,DD);
	%% build net vels
	if DD.switchs.netUstuff
		[map.AntiCycs.vel.net.mean,map.Cycs.vel.net.mean]=netVels(DD,map);
    end
    %% scatterstuff
    scat=scat{1}; 
	%% save
    save([DD.path.analyzed.name,'scat.mat'],'scat');
	save([DD.path.analyzed.name,'maps.mat'],'-struct','map');
	save([DD.path.analyzed.name,'vecs.mat'],'-struct','vecs');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinMax=resortTracks(DD,MinMax,TT)
	subfields=DD.FieldKeys.trackPlots;
	track= TT.eddy.trck;
	TT.lat=extractdeepfield(track,'geo.lat');
	TT.lon=extractdeepfield(track,'geo.lon');
	for subfield=subfields'; sub=subfield{1};
		%% nicer for plots
		collapsedField=strrep(sub,'.','');
		TT.(collapsedField) =  extractdeepfield(track,sub);
		%% get statistics for track
		[TT,MinMax]=getStats(TT,MinMax,collapsedField);
	end
	switch TT.sense
		case -1
			outfile=[DD.path.analyzedTracks.AC.name,TT.fname];
		case 1
			outfile=[DD.path.analyzedTracks.C.name,TT.fname];
	end
	%% save
	save(outfile,'-struct','TT');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V]=getVecs(eddy,V)
	V.lat=[V.lat extractdeepfield(eddy,'trck.geo.lat')];
	V.age=[V.age eddy.trck(end).age];
	death=eddy.trck(end).geo;
	V.death.lat=[V.death.lat  cat(1,death.lat)];
	V.death.lon=[ V.death.lon cat(1,death.lon)];
	birth=eddy.trck(1).geo;
	V.birth.lat=[ V.birth.lat cat(1,birth.lat)];
	V.birth.lon=[ V.birth.lon cat(1,birth.lon)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TT]=getTrack(DD,jj)
	TT.fname=DD.path.tracks.files(jj).name;
	TT.filename = [DD.path.tracks.name  TT.fname];
	try
		TT.eddy=load(TT.filename);
	catch corrupt
		warning(corrupt.identifier,corrupt.getReport)
		disp('skipping!')
		TT=[]; return
	end
	TT.sense=TT.eddy.trck(1).sense.num;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function	[strctr, eddy]=TRstructure(MAP,eddy,DD)
	switch DD.parameters.trackingRef
		case 'CenterOfVolume'
			CoV=extractdeepfield(eddy,'trck.volume.center.lin');
		case 'centroid'
			CoV=extractdeepfield(eddy,'trck.centroid.lin');
	end
	strctr.idx=MAP.idx(CoV)	;
	strctr.idxLargeMap=CoV	;
	%% delete out of bounds values
	outofbounds=isnan(strctr.idx);
	strctr.idx(outofbounds) = [];
	strctr.idxLargeMap(outofbounds) = [];
	eddy.trck(outofbounds) = [];
	tracklen=numel(eddy.trck);
	strctr.length=(1:tracklen);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [age]=TRage(MAP,eddy)
	[age,count]=protoInit(MAP.proto);
	for tt=MAP.strctr.length
		idx=MAP.strctr.idx(tt);
		count(idx)=count(idx) + 1;
		age_now=eddy.trck(tt).age;
		age.mean(idx)=meanOnFly(age.mean(idx), age_now, count(idx));
		age.std(idx)=stdOnFly(age.std(idx), age_now, count(idx));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function	amp=TRamp(MAP,eddy)
	[amp.to_mean,count]=protoInit(MAP.proto);
	[amp.to_contour,~]=protoInit(MAP.proto);
	[amp.to_ellipse,~]=protoInit(MAP.proto);
	for tt=MAP.strctr.length
		idx=MAP.strctr.idx(tt);
		count(idx)=count(idx) + 1;
		amp.to_mean.mean(idx)=meanOnFly(	amp.to_mean.mean(idx), eddy.trck(tt).peak.amp.to_mean,	count(idx));
		amp.to_contour.mean(idx)=meanOnFly(	amp.to_contour.mean(idx), eddy.trck(tt).peak.amp.to_contour,count(idx));
		amp.to_ellipse.mean(idx)=meanOnFly(	amp.to_ellipse.mean(idx), eddy.trck(tt).peak.amp.to_ellipse,count(idx));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function	radius=TRradius(MAP,eddy)
	A={'mean';'meridional';'zonal'};
	for a=A'
		[radius.(a{1}),count]=protoInit(MAP.proto);
	end
	for tt=MAP.strctr.length
		idx=MAP.strctr.idx(tt);
		count(idx)=count(idx) + 1;
		for a=A'
			radius_now=eddy.trck(tt).radius.(a{1});
			radius.(a{1}).mean(idx)=meanOnFly(radius.(a{1}).mean(idx), radius_now, count(idx));
			radius.(a{1}).std(idx)=stdOnFly(radius.(a{1}).std(idx), radius_now, count(idx));
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function	vel=TRvel(MAP,eddy)
	A={'traj';'merid';'zonal'};
	for a=A'
		[vel.(a{1}),count]=protoInit(MAP.proto);
	end
	for tt=MAP.strctr.length(1:end-1)
		idx=MAP.strctr.idx(tt);
		count(idx)=count(idx) + 1;
		for a=A'
			dist_now=eddy.dist.num.(a{1}).m(tt);
			delT= (eddy.trck(tt+1).age - eddy.trck(tt).age) * 86400;
			vel_now = dist_now/delT;
			vel.(a{1}).mean(idx)=meanOnFly(vel.(a{1}).mean(idx), vel_now, count(idx));
			vel.(a{1}).std(idx)=stdOnFly(vel.(a{1}).std(idx), vel_now, count(idx));
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function	[dist,eddy]=TRdist(MAP,eddy)
	%% set up
	[dist,count]=setupDist(MAP);
	%% calc distances
	[eddy.dist.num,eddy.dist.drct]=diststuff(field2mat(eddy.trck,'geo')');
	for tt=MAP.strctr.length
		idx=MAP.strctr.idx(tt);
		count(idx)=count(idx) + 1;
		%% traj from birth
		newValue=eddy.dist.num.traj.fromBirth(tt);
		dist.traj.fromBirth.mean(idx)=meanOnFly(dist.traj.fromBirth.mean(idx),newValue , count(idx));
		dist.traj.fromBirth.std(idx)=stdOnFly(dist.traj.fromBirth.std(idx),newValue , count(idx));
		%% traj till death
		newValue=eddy.dist.num.traj.tillDeath(tt);
		dist.traj.tillDeath.mean(idx)=meanOnFly(dist.traj.tillDeath.mean(idx),newValue , count(idx));
		dist.traj.tillDeath.std(idx)=stdOnFly(dist.traj.tillDeath.std(idx),newValue , count(idx));
		%% zonal from birth
		newValue=eddy.dist.num.zonal.fromBirth(tt);
		dist.zonal.fromBirth.mean(idx)=meanOnFly(dist.zonal.fromBirth.mean(idx),newValue , count(idx));
		dist.zonal.fromBirth.std(idx)=stdOnFly(dist.zonal.fromBirth.std(idx),newValue , count(idx));
		%% zonal till death
		newValue=eddy.dist.num.zonal.tillDeath(tt);
		dist.zonal.tillDeath.mean(idx)=meanOnFly(dist.zonal.tillDeath.mean(idx),newValue , count(idx));
		dist.zonal.tillDeath.std(idx)=stdOnFly(dist.zonal.tillDeath.std(idx),newValue , count(idx));
		%% meridional from birth
		newValue=eddy.dist.num.merid.fromBirth(tt);
		dist.merid.fromBirth.mean(idx)=meanOnFly(dist.merid.fromBirth.mean(idx),newValue , count(idx));
		dist.merid.fromBirth.std(idx)=stdOnFly(dist.merid.fromBirth.std(idx),newValue , count(idx));
		%% meridional till death
		newValue=eddy.dist.num.merid.tillDeath(tt);
		dist.merid.tillDeath.mean(idx)=meanOnFly(dist.merid.tillDeath.mean(idx),newValue , count(idx));
		dist.merid.tillDeath.std(idx)=stdOnFly(dist.merid.tillDeath.std(idx),newValue , count(idx));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [count,singlecount]=TRvisits(MAP)
	count=MAP.proto.zeros;
	singlecount=MAP.proto.zeros;
	for tt=MAP.strctr.length
		idx=MAP.strctr.idx(tt);
		count(idx)=count(idx) + 1;
	end
	sidx=unique(MAP.strctr.idx);
	singlecount(sidx)=singlecount(sidx) + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% others
function [TT,MinMax]=getStats(TT,MinMax,cf)
	%% local
	TT.max.(cf)=nanmax(TT.(cf));
	TT.min.(cf)=nanmin(TT.(cf));
	TT.median.(cf)=nanmedian(TT.(cf));
	TT.std.(cf)=nanstd(TT.(cf));
	%% global updates
	if TT.max.(cf) > MinMax.max.(cf), MinMax.max.(cf)=TT.max.(cf); end
	if TT.min.(cf) < MinMax.min.(cf), MinMax.min.(cf)=TT.min.(cf); end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MAP,V,JJ,MinMax]=initAll(DD)
	JJ=DD.threads.tracks(labindex,1):DD.threads.tracks(labindex,2);
	MAP.AntiCycs=initMAP(DD);
	MAP.Cycs=initMAP(DD);
	V.AntiCycs.age=[];V.Cycs.age=[];V.AntiCycs.lat=[];V.Cycs.lat=[];
	V.AntiCycs.birth.lat=[];V.AntiCycs.birth.lon=[];
	V.Cycs.birth.lat=[];V.Cycs.birth.lon=[];
	V.AntiCycs.death.lat=[];V.AntiCycs.death.lon=[];
	V.Cycs.death.lat=[];V.Cycs.death.lon=[];
	for subfield=DD.FieldKeys.trackPlots'; sub=subfield{1};
		collapsedField=strrep(sub,'.','');
		MinMax.max.(collapsedField)=-inf;
		MinMax.min.(collapsedField)=inf;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ro=loadRossby(DD)
	MAP=initMAP(DD);
	Ro.file=[DD.path.Rossby.name,DD.path.Rossby.files.name];
	Ro.large.radius=getfield(load([DD.path.Rossby.name 'RossbyRadius.mat']),'out');
	Ro.large.radius(Ro.large.radius==0)=nan;
	Ro.small.radius=MAP.proto.nan;
	%%
	Ro.large.phaseSpeed=getfield(load([DD.path.Rossby.name 'RossbyPhaseSpeed.mat']),'out');
	Ro.large.phaseSpeed(Ro.large.phaseSpeed==0)=nan;
	Ro.small.phaseSpeed=MAP.proto.nan;
	%% nanmean to smaller map
	lin=MAP.idx;
	T=disp_progress('init','reallocating rossby stuff indices to output map');
	uni=unique(lin(lin~=0 & ~isnan(lin)));
	for li=uni
		T=disp_progress('show',T,numel(uni),10);
		Ro.small.radius(li)=nanmean(Ro.large.radius(lin==li));
		Ro.small.phaseSpeed(li)=nanmean(Ro.large.phaseSpeed(lin==li));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinMax=globalExtr(MinMax)
	%% collect high-scores for all tracks from threads
	for cf=fieldnames(MinMax.max)'; cf=cf{1};
		MinMax.max.(cf)=gop(@max, MinMax.max.(cf)) ;
		MinMax.min.(cf)=gop(@min, MinMax.min.(cf)) ;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zonMean=zonmeans(M,DD)
	zonMean=M;
	for sense=DD.FieldKeys.senses'; sen=sense{1};
		N=M.(sen).visits.all;
		for Field=DD.FieldKeys.MeanStdFields'; field=Field{1};
			wms=weightedZonMean(cell2mat(extractdeepfield(M.(sen),field)),N);
			fields = textscan(field,'%s','Delimiter','.');
			zonMean.(sen)=setfield(zonMean.(sen),fields{1}{:},wms);
		end
	end
	if DD.switchs.RossbyStuff
		%%
		zonMean.Rossby.small.radius=nanmean(M.Rossby.small.radius,2);
		zonMean.Rossby.large.radius=nanmean(M.Rossby.large.radius,2);
		%%
		zonMean.Rossby.small.phaseSpeed=nanmean(M.Rossby.small.phaseSpeed,2);
		zonMean.Rossby.large.phaseSpeed=nanmean(M.Rossby.large.phaseSpeed,2);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT=weightedZonMean(MS,weight)
	warning('off','MATLAB:divideByZero') %#ok<*RMWRN>
	weight(weight==0)=nan;
	OUT.mean=nansum(MS.mean.*weight,2)./nansum(weight,2);
	OUT.std=nansum(MS.std.*weight,2)./nansum(weight,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,count]=setupDist(MAP)
	A={'traj';'merid';'zonal'};
	B={'fromBirth';'tillDeath'};
	for a=A'
		for b=B'
			[dist.(a{1}).(b{1}),count]=protoInit(MAP.proto);
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,drct]=diststuff(geo)
	geo=[geo(1,:); geo];
	%%
	[d.traj.deg, drct.traj]=distance(geo(1:end-1,:),geo(2:end,:));
	d.traj.m=deg2rad(d.traj.deg)*earthRadius;
	d.traj.fromBirth = cumsum(d.traj.m);
	d.traj.tillDeath = flipud(cumsum(flipud(d.traj.m)));
	%%
	latmean=mean(geo(:,1));
	[d.zonal.deg, drct.zonal]=distance(latmean,geo(1:end-1,2),latmean,geo(2:end,2));
	drct.zonal(drct.zonal<=180 & drct.zonal >= 0) = 1;
	drct.zonal(drct.zonal> 180 & drct.zonal <= 360) = -1;
	d.zonal.m=deg2rad(d.zonal.deg).*drct.zonal * earthRadius;
	d.zonal.fromBirth = cumsum(d.zonal.m);
	d.zonal.tillDeath = flipud(cumsum(flipud(d.zonal.m)));
	%%
	lonmean=mean(geo(:,2));
	[d.merid.deg, drct.merid]=distance(geo(1:end-1,1),lonmean,geo(2:end,1),lonmean);
	drct.merid(drct.merid<=90 & drct.merid >= 270) = 1;
	drct.merid (drct.merid > 90 & drct.merid < 270) = -1;
	d.merid.m=deg2rad(d.merid.deg).*drct.merid * earthRadius;
	d.merid.fromBirth = cumsum(d.merid.m);
	d.merid.tillDeath = flipud(cumsum(flipud(d.merid.m)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [param,count]=protoInit(proto,type)
	if nargin < 2, type='nan'; end
	param.mean=proto.(type);
	param.std=proto.(type);
	count=proto.zeros;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ALL=mergeMapData(MAP,DD)
	if DD.threads.num>1 && ~DD.debugmode
		ALL=spmdCase(MAP,DD);
	else
		ALL=MAP{1};
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function	vecs=mergeVecData(vecs)
	vecs=vecs{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function old=comboMS(old,new,DD)
	subfieldstrings=DD.FieldKeys.MeanStdFields;
	for ff=1:numel(subfieldstrings)
		%%	 extract current field to mean/std level
		value.new=cell2mat(extractdeepfield(new,[subfieldstrings{ff}]));
		value.old=cell2mat(extractdeepfield(old,[subfieldstrings{ff}]));
		%% nan2zero
		value.new.mean(isnan(value.new.mean))=0;
		value.old.mean(isnan(value.old.mean))=0;
		value.new.std(isnan(value.new.std))=0;
		value.old.std(isnan(value.old.std))=0;
		%% combo update
		combo.mean=ComboMean(new.visits.all,old.visits.all,value.new.mean,value.old.mean);
		combo.std=ComboStd(new.visits.all,old.visits.all,value.new.std,value.old.std);
		%% set to updated values
		fields = textscan(subfieldstrings{ff},'%s','Delimiter','.');
		meanfields={[fields{1};'mean']};
		stdfields={[fields{1};'std']};
		old=setfield(old,meanfields{1}{:},combo.mean);
		old=setfield(old,stdfields{1}{:},combo.std);
	end
	old.visits.all=old.visits.all + new.visits.all;
	old.visits.single=old.visits.single + new.visits.single;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ALL=spmdCase(MAP,DD)
	subfieldstrings=DD.FieldKeys.MeanStdFields;
	map=MAP{1};
	for sense=[{'AntiCycs'},{'Cycs'}];	sen=sense{1};
		ALL.(sen)=map.(sen);
		T=disp_progress('init',['combining results from all threads - ',sen,' ']);
		for tt=1:DD.threads.num
			T=disp_progress('calc',T,DD.threads.num,DD.threads.num);
			new = cell2mat(extractfield(MAP{tt},sen));
			if tt>1
				for ff=1:numel(subfieldstrings)
					%%	 extract current field to mean/std level
					value.new=cell2mat(extractdeepfield(new,[subfieldstrings{ff}]));
					value.old=cell2mat(extractdeepfield(old,[subfieldstrings{ff}]));
					%% nan2zero
					value.new.mean(isnan(value.new.mean))=0;
					value.old.mean(isnan(value.old.mean))=0;
					value.new.std(isnan(value.new.std))=0;
					value.old.std(isnan(value.old.std))=0;
					%% combo update
					combo.mean=ComboMean(new.visits.all,old.visits.all,value.new.mean,value.old.mean);
					combo.std=ComboStd(new.visits.all,old.visits.all,value.new.std,value.old.std);
					%% set to updated values
					fields = textscan(subfieldstrings{ff},'%s','Delimiter','.');
					meanfields={[fields{1};'mean']};
					stdfields={[fields{1};'std']};
					ALL.(sen)=setfield(ALL.(sen),meanfields{1}{:},combo.mean);
					ALL.(sen)=setfield(ALL.(sen),stdfields{1}{:},combo.std);
				end
				ALL.(sen).visits.all=ALL.(sen).visits.all + new.visits.all;
				ALL.(sen).visits.single=ALL.(sen).visits.single+ new.visits.single;
			end
			old=ALL.(sen);
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [map,vecs]=mergeThreadData(minMax,map,DD,vecs)
	try
		minMax=minMax{1};
	catch %#ok<CTCH>
		error('exit debug mode!')
	end
	map=mergeMapData(map,DD);
	vecs=mergeVecData(vecs);
	map.minMax=minMax;
	vecs.minMax=minMax;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MAP=initMAP(DD)
	MAP=load([DD.path.root,'protoMaps.mat']);
	subfieldstrings=DD.FieldKeys.MeanStdFields;
	[MSproto,~]=protoInit(MAP.proto);
	for ff=1:numel(subfieldstrings)
		fields = textscan(subfieldstrings{ff},'%s','Delimiter','.');
		meanfields={[fields{1};'mean']};
		stdfields={[fields{1};'std']};
		MAP=setfield(MAP,meanfields{1}{:},MSproto.mean);
		MAP=setfield(MAP,stdfields{1}{:},MSproto.std);
	end
	MAP.visits.single=MAP.proto.zeros;
	MAP.visits.all=MAP.proto.zeros;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
