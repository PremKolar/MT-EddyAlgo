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
	%% sequential
	seq_body(DD);
	%% update infofile
	save_info(DD);
end
%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq_body(DD)
	%% set up tracking procedure
	[tracks,OLD,phantoms]=set_up_init(DD);
	numDays=DD.checks.passed.total;
	%% start tracking
	T=disp_progress('init','tracking');
	for jj=2:numDays
		T=disp_progress('disp',T,numDays-1,499);
		%% set up current day
		[NEW]=set_up_today(DD,jj);
		%% do calculations and archivings
		[OLD,tracks]=operate_day(OLD,NEW,tracks,DD,jj,phantoms);
	end
end
function [OLD,tracks]=operate_day(OLD,NEW,tracks,DD,jj,phantoms)
	%% in case of full globe only
	if phantoms
		[NEW]=kill_phantoms(NEW);
	end
	%% find minium distances between old an new time step eddies
	[MinDists]=get_min_dists(OLD,NEW);
	%% determine which ones are tracked/died/new
	TDB=tracked_dead_born(MinDists);
	%% filter for distance per time step threshold
	TDB=filter4threshold(TDB,MinDists,DD.thresh.dist);
	%% append tracked to respective cell of temporary archive 'tracks'
	[tracks,OLD,NEW]=append_tracked(TDB,tracks,MinDists,OLD,NEW);
	%% append new ones to end of temp archive
	[tracks,NEW]=append_born(TDB, tracks, NEW);
	%% write/kill dead
	[tracks]=archive_dead(TDB, tracks, OLD.eddies, DD, jj);
	%% swap
	OLD=NEW;
	
end
%%%%%%%%%%%%%%%%%%% subs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NEW]=set_up_today(DD,jj)
	eddies=read_fields(DD,jj,'eddies');
	NEW.eddies=rmfield(eddies,'filename');
	%% get delta time
	NEW.time.daynum=DD.checks.passed.daynums(jj);
	NEW.time.delT=DD.checks.del_t(jj);
	[NEW.LON,NEW.LAT]=get_geocoor(NEW.eddies);
end
function [tracks,OLD,phantoms]=set_up_init(DD)
	%% determine whether double eddies might be present due to full lon
	phantoms=logical(cell2mat(extractdeepfield(read_fields(DD,1,'cuts'),'params.full_globe.x')));
	%% read eddies
	eddies=read_fields(DD,1,'eddies');
	[tracks,OLD.eddies]=init_day_one(eddies);
	%% append geo-coor vectors for min_dist function
	[OLD.LON,OLD.LAT]=get_geocoor(OLD.eddies);
end
function [tracks]=archive_dead(TDB, tracks, old,DD,jj)
	for sense=fieldnames(TDB)';	sen=sense{1};
		%% collect all ID's in archive
		ArchIDs=cat(2,tracks.(sen).ID);
		%% all indeces in old set of dead eddies
		dead_idxs=TDB.(sen).idx.inOld.dead;
		%% init logical of which ones to be deleted from archive
		kill_idxs=false(size(ArchIDs));
		%% loop over dead indeces
		for idx=dead_idxs;
			%% find index in archive
			AIdx = ArchIDs==old.(sen)(idx).ID;
			age = tracks.(sen)(AIdx).age;
			if age >= DD.thresh.life
				%% write to 'heap'
				archive(tracks.(sen)(AIdx).track{1}, DD.path,jj)
			end
			%% all dead get deleted
			kill_idxs(AIdx)=true;
		end
		%% kill in 'stack'
		tracks.(sen)(kill_idxs)=[];
	end
end
function archive(trck,path,jj)
	%% write out file (one per eddy)
	cc=1;
	EoD=['TRACK', sprintf('%05i',cc)];
	filename=[path.tracks.name regexprep(path.eddies.files(jj).name, 'EDDIE', EoD)];
	while true
		cc=cc+1;
		EoD=['TRACK', sprintf("%03i",cc)];
		filename=[path.tracks.name regexprep(path.eddies.files(jj).name, 'EDDIE', EoD)];
		if ~exist(filename,'file'), break; end
	end
	trck(end).filename=filename;
	save(trck(end).filename,'trck');
end
function [tracks,NEW]=append_born(TDB, tracks,NEW)
	for sense=fieldnames(TDB)';	sen=sense{1};
		maxID=max(cat(1,tracks.(sen).ID));
		NN=find(TDB.(sen).flags.inNew.born);
		if size(NN,1)>size(NN,2), NN=NN'; end
		if ~isempty(NN)
			for nn=NN
				NEW.eddies.(sen)(nn).ID=maxID+1;
				tracks.(sen)(end+1).track={NEW.eddies.(sen)(nn)}; % new position
				tracks.(sen)(end).ID=maxID+1; % IDs are unique/dont get reused
				tracks.(sen)(end).age=0;
				maxID=maxID+1;
			end
		end
	end
end
function [tracks,OLD,NEW]=append_tracked(TDB,tracks,MinDists,OLD,NEW)
	fn=fieldnames(TDB);
	for sense=fn';
		sen=cell2mat(sense);
		ArchIds=cat(2,tracks.(sen).ID);
		%% loop over successfully tracked eddies
		for nn=find(TDB.(sen).flags.inNew.tracked)
			idx=MinDists.(sen).new2old.idx(nn); % get index in old data
			ID =	OLD.eddies.(sen)(idx).ID; % get ID
			AIdx = ArchIds==ID;			% get index in archive (tracks)
			age = OLD.eddies.(sen)(idx).age; % get age
			NEW.eddies.(sen)(nn).ID= ID; % set ID accordingly for new data
			NEW.eddies.(sen)(nn).age= age + NEW.time.delT	; % set age accordingly for new data
			tracks.(sen)(AIdx).age=age + NEW.time.delT;		% update age in archive
			tracks.(sen)(AIdx).track{1}=[tracks.(sen)(AIdx).track{1}, NEW.eddies.(sen)(nn)]; % append to archive
		end
	end
end
function [tracks,new_eddies]=init_day_one(eddies)
	%% init day one
	new_eddies=rmfield(eddies,'filename');
	%% init IDs
	for sense=fieldnames(new_eddies)';	sen=sense{1};
		%% set initial ID's etc
		for ee=1:numel(new_eddies.(sen));
			new_eddies.(sen)(1,ee).ID=ee;
			tracks.(sen)(ee).track={new_eddies.(sen)(1,ee)};
			tracks.(sen)(ee).ID=ee;
			tracks.(sen)(ee).age=0;
		end
	end
end
function [TDB]=filter4threshold(TDB,MD,thresh)
	fn=fieldnames(MD);
	for sense=fn'	;	sen=cell2mat(sense);
		dist=MD.(sen).new2old.dist;
		tooQuick = (dist > thresh);
		TDB.(sen).flags.inNew.tracked(tooQuick) = false;
		TDB.(sen).flags.inNew.tooQuick = tooQuick;
		TDB.(sen).flags.inNew.born(tooQuick) = true; % all new born
	end
end
function [TDB]=tracked_dead_born(MD)
	for sense=fieldnames(MD)'	;	sen=sense{1};
		%% idx in old set of min dist claims by new set
		n2oi=MD.(sen).new2old.idx;
		%% idx in new set of min dist claims by old set
		o2ni=MD.(sen).old2new.idx;
		%% min dist values of old set of n2oi
		do=MD.(sen).old2new.dist(n2oi)';
		%% respective min dist values in new set
		dn=MD.(sen).new2old.dist;
		%% agreement among new and old ie definite tracking (with respect to new set)
		try
			TDB.(sen).flags.inNew.tracked = (do == dn);
		catch
			TDB.(sen).flags.inNew.tracked = (do == dn');
		end
		%% flag for fresh eddies with respect to new set
		TDB.(sen).flags.inNew.born = ~TDB.(sen).flags.inNew.tracked;
		%% indeces of deceised eddies with respect to old set
		TDB.(sen).idx.inOld.dead = setdiff((1:length(o2ni)),n2oi);
	end
end
function [N]=kill_phantoms(N)
	
	%step: 1/9
%0 %
%time so far:   00:00:00
%time to go:    calculating...
%Cannot remove an empty or out-of-range index from an undefined variable.

%Error in S04_track_eddies>kill_phantoms (line 201)
				%N.(field).(sen)(xx)=[];

%Error in S04_track_eddies>operate_day (line 35)
		%[NEW]=kill_phantoms(NEW);

%Error in S04_track_eddies>seq_body (line 29)
		%[OLD,tracks]=operate_day(OLD,NEW,tracks,DD,jj,phantoms);

%Error in S04_track_eddies (line 13)
	%seq_body(DD);

%Error in Sall (line 6)
%tic;S04_track_eddies;T(5).toc=toc;
 
%201 				N.(field).(sen)(xx)=[];
%K>> 

	
	
	% this takes care of potential double eddies due to map-overlap. (see S00...)
	
	fn=fieldnames(N.LON);
	for sense=fn'	;	sen=cell2mat(sense);
		%% search for identical eddies
		[LOM.a,LOM.b]=meshgrid(N.LON.(sen),N.LON.(sen));
		[LAM.a,LAM.b]=meshgrid(N.LAT.(sen),N.LAT.(sen));
		LONDIFF=abs(LOM.a - LOM.b);
		DIST=floor(real(acos(sind(LAM.a).*sind(LAM.b) + cosd(LAM.a).*cosd(LAM.b).*cosd(LONDIFF)))*earthRadius); % floor for rounding errors.. <1m -> identity
		DIST(logical(eye(size(DIST))))=nan; % nan self distance
		[Y,~]=find(DIST==0);
%% kill
		N.eddies.(sen)(Y)=[];		
		
	end
end
function [MD]=get_min_dists(OLD,NEW)
	senses=fieldnames(NEW.eddies);
	for sense=senses'	;	sen=cell2mat(sense);
		[LOM.new,LOM.old]=meshgrid(NEW.LON.(sen),OLD.LON.(sen));
		[LAM.new,LAM.old]=meshgrid(NEW.LAT.(sen),OLD.LAT.(sen));
		LONDIFF=abs(LOM.new - LOM.old);
		DIST=real(acos(sind(LAM.new).*sind(LAM.old) + cosd(LAM.new).*cosd(LAM.old).*cosd(LONDIFF)))*earthRadius;
		%% find min dists
		[MD.(sen).new2old.dist,MD.(sen).new2old.idx]=min(DIST,[],1);
		[MD.(sen).old2new.dist,MD.(sen).old2new.idx]=min(DIST,[],2);
	end
end
function [LON, LAT]=get_geocoor(eddies)
	senses=fieldnames(eddies);
	for sense=senses'
		sen=sense{1};
		LON.(sen)=extractfield(cat(1,eddies.(sen).geo),'lon');
		LAT.(sen)=extractfield(cat(1,eddies.(sen).geo),'lat');
	end
end

