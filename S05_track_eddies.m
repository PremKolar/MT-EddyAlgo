%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 17:04:31
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inter-allocate different time steps to determine tracks of eddies
function S05_track_eddies
    %% init   
    DD=initialise('eddies',mfilename);   
    %% rm old files
    rmoldtracks(DD)
    %% parallel!
    init_threads(2);
    main(DD)
    %% update infofile
    conclude(DD)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD)
    disp(['using all eddies from ' DD.path.eddies.name, ' !!!'])
    if DD.debugmode
        spmd_body(DD)
    else
        spmd(2)
            spmd_body(DD)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rmoldtracks(DD)
    if ~isempty(DD.path.tracks.files)
        if DD.overwrite
            system(['rm -r ' DD.path.tracks.name '*.mat']);
        else
            error('mv old tracks first')
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD)
    %% one thread do cycs, other acycs
        switch labindex
            case 2
                sen='cyclones';
            case 1
                sen='anticyclones';
        end
    %% set up tracking procedure
    [tracks,OLD,phantoms]=set_up_init(DD,sen);
    numDays=DD.checks.passedTotal;
    %% start tracking
    T=disp_progress('init',['tracking ' sen]);
    for jj=2:numDays
        T=disp_progress('disp',T,numDays-1,499);
        %% set up current day
        [NEW]=set_up_today(DD,jj,sen);
        %% do calculations and archivings
        [OLD,tracks]=operate_day(OLD,NEW,tracks,DD,phantoms,sen);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OLD,tracks]=operate_day(OLD,NEW,tracks,DD,phantoms,sen)
    %% in case of full globe only
    if phantoms
        [NEW]=kill_phantoms(NEW,sen);
    end
    %% find minium distances between old and new time step eddies
    [MinDists,passLog]=EligibleMinDistsMtrx(OLD,NEW,sen,DD);
    [tracks(:).passLog]=deal(passLog);
    %% determine which ones are tracked/died/new
    TDB=tracked_dead_born(MinDists,sen);
    %% append tracked to respective cell of temporary archive 'tracks'
    [tracks,NEW]=append_tracked(TDB,tracks,OLD,NEW,sen);
    %% append new ones to end of temp archive
    [tracks,NEW]=append_born(TDB, tracks, OLD,NEW,sen);
    %% write/kill dead
    [tracks]=archive_dead(TDB, tracks, OLD.eddies, DD, sen);
    %% swap
    OLD=NEW;
    
end
%%%%%%%%%%%%%%%%%%% subs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tracks,NEW]=append_tracked(TDB,tracks,OLD,NEW,sen)
    %% get
    ID.arch=cat(2,tracks.ID);
    flag.new=TDB.(sen).inNew.tracked;
    idx.old=TDB.(sen).inNew.n2oi(flag.new); % get index in old data
    ID.old =cat(2,OLD.eddies.(sen)(idx.old).ID); % get ID
    IDc=num2cell(ID.old);
    %% find position in archive
    [~,idx.arch] = ismember(ID.old,ID.arch);
    %%%%%%%%%%%%%%%%%%%%%%%%%%% TEMP SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%
    NEW.time.delT(isnan(NEW.time.delT))=round(nanmedian(NEW.time.delT));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    age = num2cell(cat(2,tracks(idx.arch).age) + NEW.time.delT); % get new age
    %% set
    [NEW.eddies.(sen)(flag.new).ID] = deal(IDc{:}); % set ID accordingly for new data
    [NEW.eddies.(sen)(flag.new).age] = deal(age{:}); % set age accordingly for new data
    [tracks(idx.arch).age]= deal(age{:});		% update age in archive
    %% append tracks into track cells
    idx.new=find(flag.new);
    for aa=1:length(idx.arch)
        aidx=idx.arch(aa);
        len=tracks(aidx).length+1;
        tracks(aidx).track{1}(len)=NEW.eddies.(sen)(idx.new(aa));
        tracks(aidx).length=len;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NEW]=set_up_today(DD,jj,sen)
    eddies=read_fields(DD,jj,'eddies');
    NEW.eddies=rmfield(eddies,'filename');
    %% get delta time
    NEW.time.daynum=DD.checks.passed(jj).daynums;
    NEW.time.delT=DD.checks.del_t(jj);
    [NEW.lon,NEW.lat]=get_geocoor(NEW.eddies,sen);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tracks,OLD,phantoms]=set_up_init(DD,sen)
    %% determine whether double eddies might be present due to full lon
    phantoms=strcmp(DD.map.window.type,'globe');
    %% read eddies
    eddies=read_fields(DD,1,'eddies');
    [tracks,OLD.eddies]=init_day_one(eddies,sen);
    %% append geo-coor vectors for min_dist function
    [OLD.lon,OLD.lat]=get_geocoor(OLD.eddies,sen);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tracks]=archive_dead(TDB, tracks, old,DD,sen)
    %% collect all ID's in archive
    ArchIDs=cat(2,tracks.ID);
    %% all indeces in old set of dead eddies
    dead_idxs=TDB.(sen).inOld.dead;
    %% find which ones to write and kill
    AIdxdead = find(ismember(ArchIDs',cat(1,old.(sen)(dead_idxs).ID)));
    age = cat(1,tracks(AIdxdead).age);
    id = cat(1,tracks(AIdxdead).ID);
    pass = age >= DD.thresh.life;
    %%  write to 'heap'
    if any(pass)
        lens=cat(2,tracks(AIdxdead(pass)).length);
        ll=0;
        for pa=find(pass)'; ll=ll+1;
            archive(tracks(AIdxdead(pa)).track{1}(1:lens(ll)), DD.path,id(pa));
        end
    end
    %% kill in 'stack'
    tracks(AIdxdead)=[];	% get rid of dead matter!
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function archive(trck,path,id)
    %% write out file (one per eddy)
    EoD=[sprintf('%07i',id)];
    filename=[ path.tracks.name 'TRACK' datestr(trck(1).daynum,'yyyymmdd')...
        '-' datestr(trck(end).daynum,'yyyymmdd') '_id' EoD '.mat'];
    trck(end).filename=filename;
    save(trck(end).filename,'trck');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tracks,NEW]=append_born(TDB, tracks,OLD,NEW,sen)
    maxID=max(max([cat(2,tracks.ID) NEW.eddies.(sen).ID  OLD.eddies.(sen).ID]));
    flag.born.inNew=(TDB.(sen).inNew.born);
    if any(flag.born.inNew)
        %% new Ids and new indices (appended to end of tracks)
        newIds=num2cell(maxID+1:maxID+sum(flag.born.inNew));
        newendIdxs=numel(tracks)+1:numel(tracks)+sum(flag.born.inNew);
        %% deal new ids to eddies
        [NEW.eddies.(sen)(flag.born.inNew).age]=deal(0);
        [NEW.eddies.(sen)(flag.born.inNew).ID]=deal(newIds{:});
        %% deal eddies to archive and pre alloc
        idx.born.inNew=find(flag.born.inNew);nn=0;
        for tt=newendIdxs; nn=nn+1;
            tracks(tt).track{1}(1)	=NEW.eddies.(sen)(idx.born.inNew(nn));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TDB]=tracked_dead_born(MD,sen)
    %% idx in old set of min dist claims by new set
    n2oi=MD.(sen).new2old.idx;
    %% idx in new set of min dist claims by old set
    o2ni=MD.(sen).old2new.idx;
    %% index in new set of claims by old set
    io=MD.(sen).old2new.idx(n2oi)';
    %% respective index in new (from new's perspective)
    in=(1:length(n2oi));
    %% min dist values of old set of n2oi
    do=MD.(sen).old2new.dist(n2oi)';
    %% respective min dist values in new set
    dn=MD.(sen).new2old.dist;
    %% matlab sets dims arbitrarily sometimes for short vecs
    if size(do)~=size(dn), do=do'; end
    if size(io)~=size(in), io=io'; end
    %% agreement among new and old ie definite tracking (with respect to new set)  NOTE: this also takes care of nan'ed dists from nanOutOfBounds() since nan~=nan !
    TDB.(sen).inNew.tracked = ((do == dn) & (io == in));
    %% flag for fresh eddies with respect to new set
    TDB.(sen).inNew.born = ~TDB.(sen).inNew.tracked;
    %% indeces of deceised eddies with respect to old set
    TDB.(sen).inOld.dead=~ismember(1:length(o2ni),n2oi(TDB.(sen).inNew.tracked));
    %% remember cross ref
    TDB.(sen).inNew.n2oi=n2oi;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inout]=kill_phantoms(inout,sen)
    %% search for identical eddies
    [LOM.a,LOM.b]=meshgrid(inout.lon.(sen),inout.lon.(sen));
    [LAM.a,LAM.b]=meshgrid(inout.lat.(sen),inout.lat.(sen));
    lonDIFF=abs(LOM.a - LOM.b);
    DIST=floor(real(acos(sind(LAM.a).*sind(LAM.b) + cosd(LAM.a).*cosd(LAM.b).*cosd(lonDIFF)))*earthRadius); % floor for rounding errors.. <1m -> identity - triu so that only one of the twins gets deleted
    DIST(logical(triu(ones(size(DIST)))))=nan;% nan self distance and lower triangle (we only need one of the two for each pair)
    [Y,~] = find(DIST==0);
    %% kill
    inout.eddies.(sen)(Y)=[];
    inout.lon.(sen)(Y)=[];
    inout.lat.(sen)(Y)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closeEnough=nanOutOfBounds(NEW,OLD)
    %% get locations of new eddies
    newLin=cat(1,NEW.trackref);
    %% get possible (future) indeces for old eddies
    oldEllipIncs=rmfield(cell2mat(extractfield(OLD,'projLocsMask')),'logical');
    %% build mask. rows -> new, cols -> old
    closeEnough=false(numel(oldEllipIncs),numel(newLin));
    for ii=1:numel(oldEllipIncs)
        closeEnough(ii,:)=ismember(cat(2,newLin(:).lin),oldEllipIncs(ii).lin');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pass=checkAmpAreaBounds(OLD,NEW,sen,ampArea)
    %% get amp and area
    amp.old=extractdeepfield(OLD.eddies.(sen),'peak.amp.to_mean');
    amp.new=extractdeepfield(NEW.eddies.(sen),'peak.amp.to_mean');
    area.old=extractdeepfield(OLD.eddies.(sen),'area.intrp');
    area.new=extractdeepfield(NEW.eddies.(sen),'area.intrp');
    %% get factors between all new and all old
    [AMP.old,AMP.new]=ndgrid(amp.old,amp.new);
    [AREA.old,AREA.new]=ndgrid(area.old,area.new);
    AMPfac=AMP.old./AMP.new;
    AREAfac=AREA.old./AREA.new;
    %% check for thresholds
    pass=(AMPfac <= ampArea(2)) & (AMPfac >= ampArea(1))...
        & (AREAfac <= ampArea(2)) & (AREAfac >= ampArea(1));
    
    %     prcnt.amp=sum(sum(AMPfac <= ampArea(2) & AMPfac >= ampArea(1)))/numel(AMPfac)*100;
    %      prcnt.area=sum(sum(AREAfac <= ampArea(2) & AREAfac >= ampArea(1)))/numel(AMPfac)*100;
    % x=1:numel(AREAfac);
    % y=[ones(1,numel(AMPfac))*ampArea(1);ones(1,numel(AMPfac))*ampArea(2)];
    % semilogy(x,sort(AMPfac(:)),'b',x,sort(AREAfac(:)),'red',x,y,'black');
    % legend(sprintf('amp (%2.0f %% passed)',prcnt.amp),sprintf('area (%2.0f %% passed)',prcnt.area),'thresholds');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [quo,pass]=checkDynamicIdentity(OLD,NEW,sen,thresh)
    old.peak2cont=extractdeepfield(OLD.eddies.(sen),'peak.amp.to_ellipse');
    old.dynRad=extractdeepfield(OLD.eddies.(sen),'peak.amp.to_ellipse');
    new.peak2cont=extractdeepfield(NEW.eddies.(sen),'peak.amp.to_ellipse');
    new.dynRad=extractdeepfield(NEW.eddies.(sen),'peak.amp.to_ellipse');
    [P2C.new,P2C.old]=meshgrid(new.peak2cont,old.peak2cont);
    [dR.new,dR.old]=meshgrid(new.dynRad,old.dynRad);
    quo.peak2cont=P2C.new./P2C.old;
    quo.dynRad=dR.new./dR.old;
    quo.combo=(abs(log(quo.peak2cont)) + abs(log(quo.dynRad)))/2;
    
    pass= quo.combo >= thresh(1) & quo.combo <= thresh(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LOM,LAM,passLog]=nanUnPassed(LOM,LAM,pass)
    onesonly=@(M) M==1;
    for ff=fieldnames(pass)';f=ff{1};
        passLog.(f)=sum(pass.(f)(:))./numel(pass.(f));
    end
    pass.all=reshape(struct2array(pass),[size(LAM.new),numel(fieldnames(pass))]);
    pass.combo=onesonly(mean(pass.all,3));
    LOM.old(~pass.combo)=nan;
    LOM.new(~pass.combo)=nan;
    LAM.old(~pass.combo)=nan;
    LAM.new(~pass.combo)=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MD,passLog]=EligibleMinDistsMtrx(OLD,NEW,sen,DD)
    
    %% build geo loc matrices
    [LOM.new,LOM.old]=meshgrid(NEW.lon.(sen),OLD.lon.(sen));
    [LAM.new,LAM.old]=meshgrid(NEW.lat.(sen),OLD.lat.(sen));
    %%
    if  DD.switchs.IdentityCheck
        [~,pass.idc]=checkDynamicIdentity(OLD,NEW,sen,DD.thresh.IdentityCheck);
    end
    %%
    if DD.switchs.AmpAreaCheck
        [pass.AmpArea]=checkAmpAreaBounds(OLD,NEW,sen,DD.thresh.ampArea);
    end
    %%
    if DD.switchs.distlimit
        [pass.ellipseDist]=nanOutOfBounds(NEW.eddies.(sen),OLD.eddies.(sen));
    end
    %%
    if exist('pass','var')
        [LOM,LAM,passLog]=nanUnPassed(LOM,LAM,pass);
    end
    %% calc distances between all from new to all from old
    lonDIFF=abs(LOM.new - LOM.old);
    DIST=real(acos(sind(LAM.new).*sind(LAM.old) + cosd(LAM.new).*cosd(LAM.old).*cosd(lonDIFF)))*earthRadius;
    %% find min dists
    [MD.(sen).new2old.dist,MD.(sen).new2old.idx]=min(DIST,[],1);
    [MD.(sen).old2new.dist,MD.(sen).old2new.idx]=min(DIST,[],2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lon, lat]=get_geocoor(eddies,sen)
    lon.(sen)=extractfield(cat(1,eddies.(sen).geo),'lon');
    lat.(sen)=extractfield(cat(1,eddies.(sen).geo),'lat');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

