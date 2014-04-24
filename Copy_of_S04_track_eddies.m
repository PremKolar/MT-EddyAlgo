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
	%% git
	%	auto_git
end
%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seq_body(DD)
	%% set up tracking proc
	[tracks,OLD,cut]=set_up_init(DD);
	numDays=DD.checks.passed.total;
	%% start tracking
	T=disp_progress('init','tracking');
	for jj=2:numDays
		T=disp_progress('disp',T,numDays-1,499);
		%% set up current day
		[NEW]=set_up_today(DD,jj,cut);
		%% do calculations and archivings
		[OLD,tracks]=operate_day(OLD,NEW,tracks,DD,jj);
	end
end
function [OLD,tracks]=operate_day(OLD,NEW,tracks,DD,jj)
	%% in case of full globe only
	if NEW.phantoms
		[NEW]=kill_phantoms(NEW);
	end
	%%
	[MinDists]=get_min_dists(OLD,NEW);
	%%
	TDB=tracked_dead_born(MinDists);
	%%
	TDB=filter4threshold(TDB,MinDists,DD.thresh.dist);
	%%
	[tracks,OLD,NEW]=append_tracked(TDB,tracks,MinDists,OLD,NEW);
	%%
	[tracks,NEW]=append_born(TDB, tracks, NEW);
	%%
	[tracks]=archive_dead(TDB, tracks, OLD.eddies, DD, jj);
	%% swap
	OLD=NEW;
	
end
%%%%%%%%%%%%%%%%%%% subs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NEW]=set_up_today(DD,jj,cut)
	eddies=read_fields(DD,jj,'eddies');
	NEW.eddies=rmfield(eddies,'filename');
	NEW.CoVlin=get_CoVlin(NEW.eddies);
	NEW.time.daynum=DD.checks.passed.daynums(jj);
	NEW.time.delT=DD.checks.del_t(jj);
	
	for sense=fieldnames(NEW.eddies)';	sen=sense{1};	
		for ee=1:length(NEW.eddies.(sen))
		NEW.eddies.(sen)(ee).length=1;
		end
	end
	
	NEW.eddies
	
	[NEW.LON,NEW.LAT]=get_geocoor(NEW.eddies);
	NEW.phantoms=cut.params.full_globe.x;
end
function [tracks,OLD,cut]=set_up_init(DD)
	cut=read_fields(DD,1,'cuts');
	eddies=read_fields(DD,1,'eddies');
	[tracks,OLD.eddies]=init_day_one(eddies);
	OLD.CoVlin=get_CoVlin(OLD.eddies);
	[OLD.LON,OLD.LAT]=get_geocoor(OLD.eddies);
end
function [tracks]=archive_dead(TDB, tracks, old,DD,jj)	
	for sense=fieldnames(TDB)';	sen=sense{1};
		ArchIDs=extractfield(tracks.(sen)(1,:),'ID');		
		for idx=TDB.(sen).idx.inOld.dead
			ID = old.(sen)(idx).ID;
			AIdx = ArchIDs==ID;
			age = tracks.(sen)(end,AIdx).age;
			if age >= DD.thresh.life
				%% write to 'heap'
				archive(tracks.(sen)(:,AIdx), DD.path,jj)
			end
			%% kill in 'stack'
			tracks.(sen)(:,AIdx)=[];	
			ArchIDs(AIdx)=[];
		end
	end
end
function archive(trck,path,jj)
	trck(end).filename=[path.tracks.name regexprep(path.eddies.files(jj).name, 'EDDIE', 'TRACK')];
	save(trck(end).filename,'trck');
end
function [tracks,NEW]=append_born(TDB, tracks,NEW)
	for sense=fieldnames(TDB)';	sen=sense{1};
		maxID=max(extractfield(tracks.(sen)(1,:),'ID'));
		for nn=find(TDB.(sen).flags.inNew.born)
			tracks.(sen)(1,end+1)=NEW.eddies.(sen)(nn);
			tracks.(sen)(1,end).ID=maxID+1;
			NEW.eddies.(sen)(nn).ID=maxID+1;
			maxID=maxID+1;
			
		end
	end
end

function [tracks,OLD,NEW]=append_tracked(TDB,tracks,MinDists,OLD,NEW)
	fn=fieldnames(TDB);
	for sense=fn';
		sen=cell2mat(sense);
		%% loop over successfully tracked eddies
		ArchIds=extractfield(tracks.(sen)(1,:),'ID');
		for nn=find(TDB.(sen).flags.inNew.tracked)
			idx=MinDists.(sen).new2old.idx(nn); % get index in old data
			ID =	OLD.eddies.(sen)(idx).ID; % get ID
			AIdx = ArchIds==ID;			% get index in archive (tracks)
			age = OLD.eddies.(sen)(idx).age; % get age in archive (tracks)
			NEW.eddies.(sen)(nn).ID= ID; % set ID accordingly for new data
			NEW.eddies.(sen)(nn).age= age + NEW.time.delT	; % set age accordingly for new data
			%tracklength=length(extractfield(tracks.(sen)(:,AIdx),'ID'));
			tracklength=tracks.(sen)(1,AIdx).length;
			tracks.(sen)(tracklength+1,AIdx)=NEW.eddies.(sen)(nn); % append to archive
			tracks.(sen)(1,AIdx).length=tracklength + 1;			
		end
	end
end
function [tracks,new_eddies]=init_day_one(eddies)
	%% init day one
	new_eddies=rmfield(eddies,'filename');
	%eddies=rmfield(eddies,'filename');
	%% init IDs
	for sense=fieldnames(new_eddies)';	sen=sense{1};
	%% set initial ID's
		for ee=1:numel(new_eddies.(sen));						
				new_eddies.(sen)(1,ee).length=1;	
			new_eddies.(sen)(1,ee).ID=ee;				
		end
	end
	tracks=new_eddies;
end
function [TDB]=filter4threshold(TDB,MD,thresh)
	fn=fieldnames(MD);
	for sense=fn'	;	sen=cell2mat(sense);
		dist=MD.(sen).new2old.dist;
		tooQuick = (dist(TDB.(sen).flags.inNew.tracked) > thresh);
		TDB.(sen).flags.inNew.tracked(tooQuick) = false;
		TDB.(sen).flags.inNew.tooQuick = tooQuick;
		TDB.(sen).flags.inNew.born(tooQuick) = true; % all new
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
		TDB.(sen).flags.inNew.tracked = (do == dn);
		%% flag for fresh eddies with respect to new set
		TDB.(sen).flags.inNew.born = ~TDB.(sen).flags.inNew.tracked;
		%% indeces of deceised eddies with respect to old set
		TDB.(sen).idx.inOld.dead = setdiff((1:length(o2ni)),n2oi);
	end
end
function [N]=kill_phantoms(N)
% this takes care of potential double eddies due to map-overlap. (see S00...)
	fn=fieldnames(N.LON);
	for sense=fn'	;	sen=cell2mat(sense);
		[LOM.a,LOM.b]=meshgrid(N.LON.(sen),N.LON.(sen));
		[LAM.a,LAM.b]=meshgrid(N.LAT.(sen),N.LAT.(sen));
		LONDIFF=abs(LOM.a - LOM.b);
		DIST=floor(real(acos(sind(LAM.a).*sind(LAM.b) + cosd(LAM.a).*cosd(LAM.b).*cosd(LONDIFF)))*earthRadius);
		DIST(logical(eye(size(DIST))))=nan; % nan self distance
		[Y,X]=find(DIST==0);
		while ~isempty(Y)
			xx=Y(end);
			%% delete respective double from y
			Y(X==xx)=[];
			if isempty(Y)
				break
			end
			%% delete phantom fields
			fn=fieldnames(N);
			for fields=fn'	;	field=cell2mat(fields);
				N.(field).(sen)(xx)=[];
			end
		end
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


function CoVlin=get_CoVlin(eddies)
	fn=fieldnames(eddies);
	for sense=fn'
		sen=cell2mat(sense);
		CoVtemp=cell2mat(extractfield(eddies.(sen),'volume'));
		CoVtemp=cell2mat(extractfield(CoVtemp,'center'));
		CoVlin.(sen)=cat(1,CoVtemp.lin);
	end
end
