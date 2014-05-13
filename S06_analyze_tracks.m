%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S06_analyze_tracks
	%% init
	DD=initialise;
	%%
	DD.threads.tracks=thread_distro(DD.threads.num,numel(DD.path.tracks.files));
	%%
	init_threads(DD.threads.num);
	spmd
		id=labindex;
		[map,vecs]=spmd_body(DD,id);
	end
	%% merge
	map=mergeMapData(map,DD);   %#ok<NASGU>
	vecs=mergeVecData(vecs);  %#ok<NASGU>
	%% save
	save([DD.path.analyzed.name,'maps.mat'],'-struct','map');
	save([DD.path.analyzed.name,'vecs.mat'],'-struct','vecs');
end

function MAP=initMAP(DD)
	MAP=load([DD.path.root,'protoMaps.mat']);
	subfieldstrings=DD.FieldKeys.MeanStdFields;
	[MSproto,~]=protoInit(MAP.proto);
	for ff=1:numel(subfieldstrings)
		fields = textscan(subfieldstrings{ff},'%s','Delimiter','.');
		meanfields={[fields{1};'mean']};
		stdfields={[fields{1};'std']};
		MAP=setfield(MAP,meanfields{1}{:},MSproto.mean)				;
		MAP=setfield(MAP,stdfields{1}{:},MSproto.std)				;
	end
	MAP.visits=MAP.proto.zeros;
	MAP.visitsSingleEddy=MAP.proto.zeros;
end

function resortTracks(DD,eddy,sense,fname)
	subfields=DD.FieldKeys.trackPlots;
	track=cell2mat(extractfield(eddy,'trck'));
	TT.lat=extractfield(cell2mat(extractfield(track,'geo')),'lat');
	TT.lon=extractfield(cell2mat(extractfield(track,'geo')),'lon');
	for subfield=subfields'; sub=subfield{1};
		collapsedField=strrep(sub,'.','');
		TT.(collapsedField) =  extractdeepfield(track,sub);
	end
	
	switch sense
		case -1
			outfile=[DD.path.analyzedTracks.AC.name,fname];
		case 1
			outfile=[DD.path.analyzedTracks.C.name,fname];
	end
	save(outfile,'-struct','TT');
	tracklistfile=[DD.path.analyzed.name, sprintf('%02i',labindex),'tracklist.txt'	];
	fid=fopen(tracklistfile,'a+');
	fprintf(fid,'%s\n',outfile);
	fclose(fid);
end

function [MAP,V]=spmd_body(DD,id)
	JJ=DD.threads.tracks(id,1):DD.threads.tracks(id,2);
	MAPac=initMAP(DD);
	MAPc=initMAP(DD);
	Vac.age=[];Vc.age=[];Vac.lat=[];Vc.lat=[];
	for jj=JJ;
		fname=DD.path.tracks.files(jj).name;
		filename = [DD.path.tracks.name  fname	];
		eddy=load(filename);
		sense=eddy.trck(1).sense.num;
		%%
		resortTracks(DD,eddy,sense,fname);
		%%
		switch sense
			case -1
				[MAPac,Vac]=MeanStdStuff(eddy,MAPac,Vac,DD);
			case 1
				[MAPc,Vc]=MeanStdStuff(eddy,MAPc,Vc,DD);
		end
		
	end
	
	MAP.AntiCycs=MAPac;
	MAP.Cycs=MAPc;
	V.AntiCycs=Vac;
	V.Cycs=Vc;
end

function [MAP,V]=MeanStdStuff(eddy,MAP,V,DD)
	MAP.strctr=TRstructure(MAP,eddy);
	[NEW.age]=TRage(MAP,eddy);
	[NEW.dist,eddy]=TRdist(MAP,eddy);
	NEW.vel=TRvel(MAP,eddy);
	NEW.radius=TRradius(MAP,eddy);
	NEW.amp=TRamp(MAP,eddy);
	[NEW.visits,NEW.visitsSingleEddy]=TRvisits(MAP);
	MAP=comboMS(MAP,NEW,DD);
	
	[V]=getVecs(eddy,V);
	
end

function	strctr=TRstructure(MAP,eddy)
	strctr.length=cell(numel(eddy),1);
	strctr.idx=cell(numel(eddy),1);
	strctr.lengthTotal=0;
	tracklen=numel(eddy.trck);
	strctr.lengthTotal=strctr.lengthTotal + tracklen;
	strctr.length=(1:tracklen);
	strctr.idx=nan(1,tracklen);
	for tt=1:tracklen
		strctr.idx(tt)=MAP.idx(eddy.trck(tt).volume.center.lin)	;
	end
end

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

function	amp=TRamp(MAP,eddy)
	
	[amp.to_mean.of_contour,count]=protoInit(MAP.proto);
	[amp.to_contour,~]=protoInit(MAP.proto);
	[amp.to_ellipse,~]=protoInit(MAP.proto);
	
	
	for tt=MAP.strctr.length
		idx=MAP.strctr.idx(tt);
		count(idx)=count(idx) + 1;
		amp.to_mean.of_contour.mean(idx)=meanOnFly(	amp.to_mean.of_contour.mean(idx), eddy.trck(tt).peak.amp.to_mean.of_contour,	count(idx));
		amp.to_mean.of_contour.mean(idx)=meanOnFly(	amp.to_mean.of_contour.mean(idx), eddy.trck(tt).peak.amp.to_contour,count(idx));
		amp.to_mean.of_contour.mean(idx)=meanOnFly(	amp.to_mean.of_contour.mean(idx), eddy.trck(tt).peak.amp.to_ellipse,count(idx));
	end
	
end

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

function	[dist,eddy]=TRdist(MAP,eddy)
	%% set up
	A={'traj';'merid';'zonal'};
	B={'fromBirth';'tillDeath'};
	for a=A'
		for b=B'
			[dist.(a{1}).(b{1}),count]=protoInit(MAP.proto);
		end
	end
	%%
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

function [param,count]=protoInit(proto,type)
	if nargin < 2, type='nan'; end
	param.mean=proto.(type);
	param.std=proto.(type);
	count=proto.zeros;
end

function ALL=mergeMapData(MAP,DD)
	if DD.threads.num>1
		ALL=spmdCase(MAP,DD);
	else
		ALL=MAP{1};
	end
end

function	vecs=mergeVecData(vecs)
	vecs=vecs{1};
end

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
		combo.mean=ComboMean(new.visits,old.visits,value.new.mean,value.old.mean);
		combo.std=ComboStd(new.visits,old.visits,value.new.std,value.old.std);
		%% set to updated values
		fields = textscan(subfieldstrings{ff},'%s','Delimiter','.');
		meanfields={[fields{1};'mean']};
		stdfields={[fields{1};'std']};
		old=setfield(old,meanfields{1}{:},combo.mean)				;
		old=setfield(old,stdfields{1}{:},combo.std)				;
		
	end
	old.visits=old.visits + new.visits;
	old.visitsSingleEddy=old.visitsSingleEddy + new.visitsSingleEddy;
end

function ALL=spmdCase(MAP,DD)
	subfieldstrings=DD.FieldKeys.MeanStdFields;
	map=MAP{1};
	for sense=[{'AntiCycs'},{'Cycs'}];	sen=sense{1};
		ALL.(sen)=map.(sen);
		T=disp_progress('init',['combining results from all threads - ',sen,' ']);
		for tt=1:DD.threads.num
			T=disp_progress('calc',T,DD.threads.num,DD.threads.num);
			new = MAP{tt};
			new=new.(sen);
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
					combo.mean=ComboMean(new.visits,old.visits,value.new.mean,value.old.mean);
					combo.std=ComboStd(new.visits,old.visits,value.new.std,value.old.std);
					%% set to updated values
					fields = textscan(subfieldstrings{ff},'%s','Delimiter','.');
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

