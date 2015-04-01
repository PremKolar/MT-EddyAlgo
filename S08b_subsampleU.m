function S08b_subsampleU
    DD = getfield(load('S09main.mat'),'DD');
    senses = DD.FieldKeys.senses;
    yearA = 1995;
    yearB = 2005;
    %%
    tl = {'velPP';'lat';'lon'};
    toLoad(numel(tl)).name=struct;
    [toLoad(:).name] = deal(tl{:});
    %%
    parfor yy = yearA:yearB
        doYear(yy,senses,DD,toLoad);
    end
    %%
    kk=0;
    for yy = yearA:yearB
        kk=kk+1;
        drawYears(DD,yy,kk,yearB,yearA);
    end
    
    ylabel('-u')
    xlabel('latitude')
    title(sprintf('yearly subsamples %d - %d',yearA,yearB));
    grid minor
    outfname = 'Usubsampled',
    savefig(DD.path.plots,72,400,300,outfname,'dpdf',DD2info(DD),14);
    
    
end
% -------------------------------------------------------------------------
function drawYears(DD,yy,kk,yearB,yearA)
    track = getfield(load(sprintf('%ssubVelYear%d.mat',DD.path.analyzed.name,yy)),'tout');
    latround = round(track.lat);
    laU = unique(latround);
    velDraw = nan(size(laU));
    cc=1;
    for la = laU
        todrawidx = latround == la;
        if sum(todrawidx)>100 && abs(la)>10
            velDraw(cc) = nanmean(track.vel(todrawidx));
        end
        cc=cc+1;
    end
    plot(laU,-velDraw,'color',rainbow(1,1,1,kk,yearB-yearA+1))
    hold on
    axis tight
end
% -------------------------------------------------------------------------
function doYear(yy,senses,DD,toLoad)
    fprintf('year %d',yy);
    outfname = sprintf('%ssubVelYear%d.mat',DD.path.analyzed.name,yy);
    if exist(outfname,'file')
        return
    end
    %%
    pattern = ['*TRACK*-' num2str(yy) '*_id*'];
    TRACKSa = dir2([DD.path.analyzedTracks.(senses{1}).name pattern]);
    TRACKSb = dir2([DD.path.analyzedTracks.(senses{2}).name pattern]);
    TRACKS = [TRACKSa(3:end); TRACKSb(3:end)];
    %     TRACKS = [TRACKSa(3:10); TRACKSb(3:10)];
    track(numel(TRACKS)).data = struct;
    %%
    for tt = 1:numel(TRACKS)
        tmp = load(TRACKS(tt).fullname,toLoad(:).name);
        tmp.vel = tmp.velPP.zonal.v';
        track(tt).data = rmfield(tmp,'velPP');
    end
    %%
    tmp = cat(2,track.data);
    tout.vel = cat(2,tmp.vel);
    tout.lat = cat(2,tmp.lat);
    tout.lon = cat(2,tmp.lon); %#ok<STRNU>
    %% save
    save(outfname,'tout')
end
% -------------------------------------------------------------------------