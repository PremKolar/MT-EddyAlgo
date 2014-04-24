%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S06_analyze_tracks
	%% init
	DD=initialise('tracks');
	%%
	DD.threads.lims_eddies=thread_distro(DD.threads.num,numel(DD.path.tracks.files));
	%%
	init_threads(DD.threads.num);
	spmd
	id=labindex;
	[map,tracks]=spmd_body(DD,id);
	end
	%% merge
	MAP=mergeMapData(map,DD);  %#ok<NASGU>
	TRACKS=mergeTracksData(tracks,DD); %#ok<NASGU>
	
	%% save
	save([DD.path.analyzed.name,'maps.mat'],'-struct','MAP');
	save([DD.path.analyzed.name,'tracks.mat'],'-struct','TRACKS');
	%% update infofile
	%	save_info(DD)
end

function	TRACKS=mergeTracksData(tracks,DD)
	if DD.threads.num>1
		TRACKS=tracks{1}; %already joined TrackData4Plot
	else
		TRACKS=tracks;
	end
end

function [MeanStd,tracks]=spmd_body(DD,id)
	%% init
	[AntiCycs,Cycs]=initACandC(DD,id);
	%% put tracks into better plottable struct
	[~,tracks.Cycs]=TrackData4Plot(Cycs,DD);
	[~,tracks.AntiCycs]=TrackData4Plot(AntiCycs,DD)	;
	%% Mean and STD maps
	MAP=load([DD.path.root,'protoMaps.mat']);
	MeanStd.AntiCycs=MeanStdStuff(AntiCycs,MAP);
	MeanStd.Cycs=MeanStdStuff(Cycs,MAP);
end


function [ALL,tracks]=TrackData4Plot(eddies,DD)
	disp('formating 4 plots..')
	%% get keys
	ALL=struct;
	tracks(numel(eddies))=struct;
	subfields=DD.FieldKeys.trackPlots;
	%% init
	for ee=1:numel(eddies)
		ALL.lat=extractfield(cell2mat(extractfield(eddies(1).track,'geo')),'lat');
		ALL.lon=extractfield(cell2mat(extractfield(eddies(1).track,'geo')),'lon');
		
		for subfield=subfields'; sub=subfield{1};
			collapsedField=strrep(sub,'.','');
			ALL.(collapsedField) =	extractdeepfield(eddies(1).track,sub);
		end
	end
	%% append
	for ee=2:numel(eddies)
		ALL.lat=[ALL.lat, extractfield(cell2mat(extractfield(eddies(ee).track,'geo')),'lat')];
		ALL.lon=[ALL.lon, extractfield(cell2mat(extractfield(eddies(ee).track,'geo')),'lon')];
		tracks(ee).lat=extractfield(cell2mat(extractfield(eddies(ee).track,'geo')),'lat');
		tracks(ee).lon=extractfield(cell2mat(extractfield(eddies(ee).track,'geo')),'lon');
		%%
		for subfield=subfields'; sub=subfield{1};
			collapsedField=strrep(sub,'.','');
			tracks(ee).(collapsedField) =  extractdeepfield(eddies(ee).track,sub);
			ALL.(collapsedField) =	[ALL.(collapsedField), tracks(ee).(collapsedField)];
		end
	end
	
	%%
	disp('sending to master..')
	ALL=gcat(ALL,2,1);
	tracks=gcat(tracks,2,1);
end


function [AntiCycs,Cycs]=initACandC(DD,id)
	JJ=DD.threads.lims_eddies(id,1):DD.threads.lims_eddies(id,2);
	Cycs(numel(JJ),1)=struct; %pre allocate
	AntiCycs(numel(JJ),1)=struct;
	ac=0; cc=0;
	T=disp_progress('init','collecting eddies');
	%%
	for jj=JJ;
		T=disp_progress('calc',T,numel(JJ),10);
		filename = [DD.path.tracks.name  DD.path.tracks.files(jj).name	];
		eddy=load(filename);
		sense=eddy.trck(1).sense.num;
		switch sense
			case 1
				ac=ac+1;
				AntiCycs(ac).track=cell2mat(extractfield(eddy,'trck'));
			case -1
				cc=cc+1;
				Cycs(cc).track=cell2mat(extractfield(eddy,'trck'));
		end
	end
	%% kill empty pre allocated space
	AntiCycs(ac+1:end)=[];
	Cycs(cc+1:end)=[];
end
function MAP=MeanStdStuff(eddies,MAP)
	MAP.strctr=TRstructure(MAP,eddies);
	disp('counting visits')
	[MAP.visits,MAP.visitsSingleEddy]=TRvisits(MAP,eddies);
	disp('age stuff')
	MAP.age=TRage(MAP,eddies);
	disp('sense stuff')
	MAP.sense=TRsense(MAP,eddies);
	disp('distance stuff')
	[MAP.dist,eddies]=TRdist(MAP,eddies);
	disp('velocity stuff')
	MAP.vel=TRvel(MAP,eddies);
	disp('radius stuff')
	MAP.radius=TRradius(MAP,eddies);
	disp('amp stuff')
	MAP.amp=TRamp(MAP,eddies);
	
	
end

function	amp=TRamp(MAP,eddies)
	
	[amp.to_mean.of_contour,count]=protoInit(MAP.proto);
	[amp.to_contour,~]=protoInit(MAP.proto);
	[amp.to_ellipse,~]=protoInit(MAP.proto);
	
	for ee=1:numel(eddies)
		for tt=MAP.strctr.length{ee}
			idx=MAP.strctr.idx{ee}(tt);
			count(idx)=count(idx) + 1;
			amp.to_mean.of_contour.mean(idx)=meanOnFly(	amp.to_mean.of_contour.mean(idx), eddies(ee).track(tt).peak.amp.to_mean.of_contour,	count(idx));
			amp.to_mean.of_contour.mean(idx)=meanOnFly(	amp.to_mean.of_contour.mean(idx), eddies(ee).track(tt).peak.amp.to_contour,count(idx));
			amp.to_mean.of_contour.mean(idx)=meanOnFly(	amp.to_mean.of_contour.mean(idx), eddies(ee).track(tt).peak.amp.to_ellipse,count(idx));
		end
	end
end


function	radius=TRradius(MAP,eddies)
	A={'mean';'meridional';'zonal'};
	for a=A'
		[radius.(a{1}),count]=protoInit(MAP.proto);
	end
	for ee=1:numel(eddies)
		for tt=MAP.strctr.length{ee}
			idx=MAP.strctr.idx{ee}(tt);
			count(idx)=count(idx) + 1;
			for a=A'
				radius_now=eddies(ee).track(tt).radius.(a{1});
				radius.(a{1}).mean(idx)=meanOnFly(radius.(a{1}).mean(idx), radius_now, count(idx));
				radius.(a{1}).std(idx)=stdOnFly(radius.(a{1}).std(idx), radius_now, count(idx));
			end
		end
	end
end
function	vel=TRvel(MAP,eddies)
	A={'traj';'merid';'zonal'};
	for a=A'
		[vel.(a{1}),count]=protoInit(MAP.proto);
	end
	for ee=1:numel(eddies)
		for tt=MAP.strctr.length{ee}(1:end-1)
			idx=MAP.strctr.idx{ee}(tt);
			count(idx)=count(idx) + 1;
			for a=A'
				dist_now=eddies(ee).dist.num.(a{1}).m(tt);
				delT= (eddies(ee).track(tt+1).age - eddies(ee).track(tt).age) * 86400;
				vel_now = dist_now/delT;
				vel.(a{1}).mean(idx)=meanOnFly(vel.(a{1}).mean(idx), vel_now, count(idx));
				vel.(a{1}).std(idx)=stdOnFly(vel.(a{1}).std(idx), vel_now, count(idx));
			end
		end
	end
end
function	[dist,eddies]=TRdist(MAP,eddies)
	%% set up
	A={'traj';'merid';'zonal'};
	B={'fromBirth';'tillDeath'};
	for a=A'
		for b=B'
			[dist.(a{1}).(b{1}),count]=protoInit(MAP.proto);
		end
	end
	%%
	for ee=1:numel(eddies)
		%% calc distances
		
		[eddies(ee).dist.num,eddies(ee).dist.drct]=diststuff(field2mat(eddies(ee).track,'geo')');
		for tt=MAP.strctr.length{ee}(1:end-1)
			idx=MAP.strctr.idx{ee}(tt);
			count(idx)=count(idx) + 1;
			%% traj from birth
			newValue=eddies(ee).dist.num.traj.fromBirth(tt);
			dist.traj.fromBirth.mean(idx)=meanOnFly(dist.traj.fromBirth.mean(idx),newValue , count(idx));
			dist.traj.fromBirth.std(idx)=stdOnFly(dist.traj.fromBirth.std(idx),newValue , count(idx));
			%% traj till death
			newValue=eddies(ee).dist.num.traj.tillDeath(tt);
			dist.traj.tillDeath.mean(idx)=meanOnFly(dist.traj.tillDeath.mean(idx),newValue , count(idx));
			dist.traj.tillDeath.std(idx)=stdOnFly(dist.traj.tillDeath.std(idx),newValue , count(idx));
			%% zonal from birth
			newValue=eddies(ee).dist.num.zonal.fromBirth(tt);
			dist.zonal.fromBirth.mean(idx)=meanOnFly(dist.zonal.fromBirth.mean(idx),newValue , count(idx));
			dist.zonal.fromBirth.std(idx)=stdOnFly(dist.zonal.fromBirth.std(idx),newValue , count(idx));
			%% zonal till death
			newValue=eddies(ee).dist.num.zonal.tillDeath(tt);
			dist.zonal.tillDeath.mean(idx)=meanOnFly(dist.zonal.tillDeath.mean(idx),newValue , count(idx));
			dist.zonal.tillDeath.std(idx)=stdOnFly(dist.zonal.tillDeath.std(idx),newValue , count(idx));
			%% meridional from birth
			newValue=eddies(ee).dist.num.merid.fromBirth(tt);
			dist.merid.fromBirth.mean(idx)=meanOnFly(dist.merid.fromBirth.mean(idx),newValue , count(idx));
			dist.merid.fromBirth.std(idx)=stdOnFly(dist.merid.fromBirth.std(idx),newValue , count(idx));
			%% meridional till death
			newValue=eddies(ee).dist.num.merid.tillDeath(tt);
			dist.merid.tillDeath.mean(idx)=meanOnFly(dist.merid.tillDeath.mean(idx),newValue , count(idx));
			dist.merid.tillDeath.std(idx)=stdOnFly(dist.merid.tillDeath.std(idx),newValue , count(idx));
		end
	end
end
function [d,drct]=diststuff(geo)
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
function	sense=TRsense(MAP,eddies)
	[sense,count]=protoInit(MAP.proto);
	for ee=1:numel(eddies)
		for tt=MAP.strctr.length{ee}
			idx=MAP.strctr.idx{ee}(tt);
			count(idx)=count(idx) + 1;
			sense_now=eddies(ee).track(tt).sense.num;
			sense.mean(idx)=meanOnFly(sense.mean(idx), sense_now, count(idx));
			sense.std(idx)=stdOnFly(sense.std(idx), sense_now, count(idx));
		end
	end
end
function [count,singlecount]=TRvisits(MAP,eddies)
	count=MAP.proto.zeros;
	singlecount=MAP.proto.zeros;
	for ee=1:numel(eddies)
		for tt=MAP.strctr.length{ee}
			idx=MAP.strctr.idx{ee}(tt);
			count(idx)=count(idx) + 1;
		end
		sidx=unique(MAP.strctr.idx{ee});
		singlecount(sidx)=singlecount(sidx) + 1;
	end
end
function age=TRage(MAP,eddies)
	[age,count]=protoInit(MAP.proto);
	for ee=1:numel(eddies)
		for tt=MAP.strctr.length{ee}
			idx=MAP.strctr.idx{ee}(tt);
			count(idx)=count(idx) + 1;
			age_now=eddies(ee).track(tt).age;
			age.mean(idx)=meanOnFly(age.mean(idx), age_now, count(idx));
			age.std(idx)=stdOnFly(age.std(idx), age_now, count(idx));
		end
	end
end
function [param,count]=protoInit(proto,type)
	if nargin < 2, type='nan'; end
	param.mean=proto.(type);
	param.std=proto.(type);
	count=proto.zeros;
end
function	strctr=TRstructure(MAP,eddies)
	strctr.length=cell(numel(eddies),1);
	strctr.idx=cell(numel(eddies),1);
	for ee=1:numel(eddies)
		tracklen=numel(eddies(ee).track);
		strctr.length{ee}=(1:tracklen);
		strctr.idx{ee}=nan(1,tracklen);
		for tt=1:tracklen
			strctr.idx{ee}(tt)=MAP.idx(eddies(ee).track(tt).volume.center.lin)	;
		end
	end
end
function ALL=mergeMapData(MAP,DD)
	if DD.threads.num>1
		ALL=spmdCase(MAP,DD);
	else
		ALL=MAP;
	end
end
function ALL=spmdCase(MAP,DD)
	subfieldstrings=DD.FieldKeys.MeanStdFields;
	T=disp_progress('init','combining results from all threads');
	map=MAP{1};
	for sense=[{'AntiCycs'},{'Cycs'}];	sen=sense{1};
		ALL.(sen)=map.(sen);
		for tt=1:DD.threads.num
			T=disp_progress('calc',T,DD.threads.num,DD.threads.num);
			new = MAP{tt};
			new=new.(sen);
			if tt>1
				for ff=1:numel(subfieldstrings)
					%%
					[value.new.mean,fields]=getSubField([subfieldstrings{ff},'.mean'],new);
					value.old.mean=getSubField([subfieldstrings{ff},'.mean'],old);
					value.new.std=getSubField([subfieldstrings{ff},'.std'],new);
					value.old.std=getSubField([subfieldstrings{ff},'.std'],old);
					%%
					value.new.mean(isnan(value.new.mean))=0;
					value.old.mean(isnan(value.old.mean))=0;
					value.new.std(isnan(value.new.std))=0;
					value.old.std(isnan(value.old.std))=0;
					%%
					combo.mean=ComboMean(new.visits,old.visits,value.new.mean,value.old.mean);
					combo.std=ComboStd(new.visits,old.visits,value.new.std,value.old.std);
					meanfields={['mean';fields{1}]};
					stdfields={['std';fields{1}]};
					ALL.(sen)=setfield(ALL.(sen),meanfields{1}{:},combo.mean)				;
					ALL.(sen)=setfield(ALL.(sen),stdfields{1}{:},combo.std)				;
				end
				ALL.(sen).visits=ALL.(sen).visits + new.visits;
			end
			old=ALL.(sen);
		end
	end
end





