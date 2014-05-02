%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 17:04:31
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inter-allocate different time steps to determine tracks of eddies
% entirely sequential no spmd...
function S04_track_eddies
	%% init
	DD=initialise('eddies');
	%% parallel!
	init_threads(2);
	spmd(2);
		seq_body(DD);
	end
	%% update infofile
	save_info(DD);
end
%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq_body(DD)
	%% one thread do cycs, other acycs
	switch labindex
		case 1
			sen='cyclones';
		case 2
			sen='anticyclones';
	end
	%% set up tracking procedure
	[tracks,OLD,phantoms]=set_up_init(DD,sen);
	numDays=DD.checks.passed.total;
	%% start tracking
	T=disp_progress('init','tracking');
	for jj=2:numDays
		T=disp_progress('disp',T,numDays-1,499);
		%% set up current day
		[NEW]=set_up_today(DD,jj,sen);
		%% do calculations and archivings
		[OLD,tracks]=operate_day(OLD,NEW,tracks,DD,jj,phantoms,sen);
	end
end
function [OLD,tracks]=operate_day(OLD,NEW,tracks,DD,jj,phantoms,sen)
	%% in case of full globe only
	if phantoms
		[NEW]=kill_phantoms(NEW,sen);
	end
	%% find minium distances between old and new time step eddies
	[MinDists]=get_min_dists(OLD,NEW,sen);
	%% determine which ones are tracked/died/new
	TDB=tracked_dead_born(MinDists,sen);
	%% filter for distance per day threshold
	dist_thresh=DD.checks.del_t(jj)*DD.thresh.dist;
	TDB=filter4threshold(TDB,MinDists,dist_thresh,sen);
	%% append tracked to respective cell of temporary archive 'tracks'
	[tracks,NEW]=append_tracked(TDB,tracks,OLD,NEW,sen);
	%% append new ones to end of temp archive
	[tracks,NEW]=append_born(TDB, tracks, NEW,sen);
	%% write/kill dead
	[tracks]=archive_dead(TDB, tracks, OLD.eddies, DD, jj,sen);
	%% swap
	OLD=NEW;
	
end
%%%%%%%%%%%%%%%%%%% subs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NEW]=set_up_today(DD,jj,sen)
	eddies=read_fields(DD,jj,'eddies');
	NEW.eddies=rmfield(eddies,'filename');
	%% get delta time
	NEW.time.daynum=DD.checks.passed.daynums(jj);
	NEW.time.delT=DD.checks.del_t(jj);
	[NEW.LON,NEW.LAT]=get_geocoor(NEW.eddies,sen);
end
function [tracks,OLD,phantoms]=set_up_init(DD,sen)
	%% determine whether double eddies might be present due to full lon
	phantoms=logical(cell2mat(extractdeepfield(read_fields(DD,1,'cuts'),'params.full_globe.x')));
	%% read eddies
	eddies=read_fields(DD,1,'eddies');
	[tracks,OLD.eddies]=init_day_one(eddies,sen);
	%% append geo-coor vectors for min_dist function
	[OLD.LON,OLD.LAT]=get_geocoor(OLD.eddies,sen);
end

function [tracks]=archive_dead(TDB, tracks, old,DD,jj,sen)
	
	%% collect all ID's in archive
	ArchIDs=cat(2,tracks.ID);
	%% all indeces in old set of dead eddies
	dead_idxs=TDB.(sen).inOld.dead;
	%% find which ones to write and kill
	AIdxdead = find(ismember(ArchIDs,cat(1,old.(sen)(dead_idxs).ID)));
	age = cat(1,tracks(AIdxdead).age);
	id = cat(1,tracks(AIdxdead).ID);
	pass = age >= DD.thresh.life;
	%%  write to 'heap'
	if any(pass)
		lens=cat(2,tracks(AIdxdead(pass)).length);
		ll=0;
		for pa=find(pass)'; ll=ll+1;
			
			archive(tracks(AIdxdead(pa)).track{1}(1:lens(ll)), DD.path,jj,id(pa));
		end
	end
	%% kill in 'stack'
	tracks(AIdxdead)=[];	% get rid of dead matter!
	
end
function archive(trck,path,jj,id)
	%% write out file (one per eddy)
	EoD=['TRACK', sprintf('%06i',id)];
	filename=[ path.tracks.name EoD regexprep(path.eddies.files(jj).name, 'EDDIE', '')];
	trck(end).filename=filename;
	save(trck(end).filename,'trck');
end
function [tracks,NEW]=append_born(TDB, tracks,NEW,sen)
	maxID=max(max([cat(2,tracks.ID) NEW.eddies.(sen).ID]));
	NN=(TDB.(sen).inNew.born);
	if any(NN)
		%% new Ids and new indices (appended to end of tracks)
		newIds=num2cell(maxID+1:maxID+sum(NN));
		newendIdxs=numel(tracks)+1:numel(tracks)+sum(NN);
		%% deal new ids to eddies
		[NEW.eddies.(sen)(NN).age]=deal(0);
		[NEW.eddies.(sen)(NN).ID]=deal(newIds{:});
		%% deal eddies to archive and pre alloc
		NN=find(NN);nn=0;
		for tt=newendIdxs; nn=nn+1;
			tracks(tt).track{1}(1)	=NEW.eddies.(sen)(NN(nn));
			tracks(tt).track{1}(30)	=tracks(tt).track{1}(1);
		end
		%% set all ages 0
		[tracks(newendIdxs).age]=deal(0);
		%% deal new ids to tracks
		[tracks(newendIdxs).ID]=deal(newIds{:});
		%% init length
		[tracks(newendIdxs).length]=deal(1);
		
	end
end
function [tracks,NEW]=append_tracked(TDB,tracks,OLD,NEW,sen)
	%% get
	ArchIds=cat(2,tracks.ID);
	NN=TDB.(sen).inNew.tracked;
	idx=TDB.(sen).inNew.n2oi(NN); % get index in old data
	ID =	cat(2,OLD.eddies.(sen)(idx).ID); % get ID
	IDc=num2cell(ID);
	%% find position in archive	
	[~,AIdx] = ismember(ID,ArchIds);
	%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMP SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%
	if any(AIdx==0)
		wrong=find(AIdx);
		NN(wrong)=[];
		IDc(wrong)=[];
		AIdx(wrong)=[];
		%	error('something wrong went wrong!!!');
		warning('something wrong went wrong!!!'); %#ok<WNTAG>
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMP SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%
	NEW.time.delT(isnan(NEW.time.delT))=1;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	age = num2cell(cat(2,tracks(AIdx).age) + NEW.time.delT); % get new age
	%% set
	[NEW.eddies.(sen)(NN).ID] = deal(IDc{:}); % set ID accordingly for new data
	[NEW.eddies.(sen)(NN).age] = deal(age{:}); % set age accordingly for new data
	[tracks(AIdx).age]= deal(age{:});		% update age in archive
	%% append tracks into track cells
	NNf=find(NN);
	for aa=1:length(AIdx)
		aidx=AIdx(aa);
		len=tracks(aidx).length+1;
		tracks(aidx).track{1}(len)=NEW.eddies.(sen)(NNf(aa));
		tracks(aidx).length=len;
	end
end
function [tracks,new_eddies]=init_day_one(eddies,sen)
	%% init day one
	new_eddies=rmfield(eddies,'filename');
	%% set initial ID's etc
	ee=(1:numel(new_eddies.(sen)));
	eec=num2cell(ee);
	[new_eddies.(sen)(1,ee).ID]=deal(eec{:});
	%% store tracks in cells (to allow for arbitr. lengths)
	edsArray=arrayfun(@(x) ({{x}}),new_eddies.(sen)(1,ee));
	[tracks(ee).track]=deal(edsArray{:});
	[tracks(ee).ID]=deal(eec{:});
	[tracks(ee).age]=deal(0);
	[tracks(ee).length]=deal(1);
	
	%% prealloc for speed
	for tt=1:length(tracks)
		tracks(tt).track{1}(30)	=tracks(tt).track{1}(1);
	end
	
end
function [TDB]=filter4threshold(TDB,MD,thresh,sen)
	dist=MD.(sen).new2old.dist;
	%% find those that were supposedly tracked, yet dont fullfill threshold (for new)
	tooQuick = ((dist > thresh) & TDB.(sen).inNew.tracked);
	%% set them to ~tracked
	TDB.(sen).inNew.tracked(tooQuick) = false;
	%% instead set them to 'born'
	TDB.(sen).inNew.born(tooQuick) = true; % all new born
	%% add respective indices for old set to 'dead' flags (track broke off)
	TDB.(sen).inOld.dead(TDB.(sen).inNew.n2oi(tooQuick))=true; % not tracked!
end
function [TDB]=tracked_dead_born(MD,sen)
	%% idx in old set of min dist claims by new set
	n2oi=MD.(sen).new2old.idx;
	%% idx in new set of min dist claims by old set
	o2ni=MD.(sen).old2new.idx;
	%% min dist values of old set of n2oi
	do=MD.(sen).old2new.dist(n2oi)';
	%% respective min dist values in new set
	dn=MD.(sen).new2old.dist;
	%% agreement among new and old ie definite tracking (with respect to new set)
	TDB.(sen).inNew.tracked = (do == dn);
	%% flag for fresh eddies with respect to new set
	TDB.(sen).inNew.born = ~TDB.(sen).inNew.tracked;
	%% indeces of deceised eddies with respect to old set
	TDB.(sen).inOld.dead=~ismember(1:length(o2ni),n2oi(TDB.(sen).inNew.tracked));
	%% remember cross ref
	TDB.(sen).inNew.n2oi=n2oi;
end
function [N]=kill_phantoms(N,sen)
	%% search for identical eddies
	[LOM.a,LOM.b]=meshgrid(N.LON.(sen),N.LON.(sen));
	[LAM.a,LAM.b]=meshgrid(N.LAT.(sen),N.LAT.(sen));
	LONDIFF=abs(LOM.a - LOM.b);
	DIST=floor(real(acos(sind(LAM.a).*sind(LAM.b) + cosd(LAM.a).*cosd(LAM.b).*cosd(LONDIFF)))*earthRadius); % floor for rounding errors.. <1m -> identity
	DIST(logical(eye(size(DIST))))=nan; % nan self distance
	Y = DIST==0;
	%% kill
	N.eddies.(sen)(Y)=[];
	N.time.(sen)(Y)=[];
	N.LON.(sen)(Y)=[];
	N.LAT.(sen)(Y)=[];
end
function [MD]=get_min_dists(OLD,NEW,sen)
	[LOM.new,LOM.old]=meshgrid(NEW.LON.(sen),OLD.LON.(sen));
	[LAM.new,LAM.old]=meshgrid(NEW.LAT.(sen),OLD.LAT.(sen));
	LONDIFF=abs(LOM.new - LOM.old);
	DIST=real(acos(sind(LAM.new).*sind(LAM.old) + cosd(LAM.new).*cosd(LAM.old).*cosd(LONDIFF)))*earthRadius;
	%% find min dists
	[MD.(sen).new2old.dist,MD.(sen).new2old.idx]=min(DIST,[],1);
	[MD.(sen).old2new.dist,MD.(sen).old2new.idx]=min(DIST,[],2);
end
function [LON, LAT]=get_geocoor(eddies,sen)
	LON.(sen)=extractfield(cat(1,eddies.(sen).geo),'lon');
	LAT.(sen)=extractfield(cat(1,eddies.(sen).geo),'lat');
end

