function [DD]=get_input
    sprintf(['\n getting user input...']);
    %%
    DD=evalUserInput;
    %%
    DD.time=catstruct(DD.time, timestuff(DD.time));
    %%
    sprintf(['\n scanning data...']);
    DD.path=catstruct(DD.path,findfiles(DD));
    %%
    sprintf(['\n setting internal parameters...']);
    [DD.pattern, DD.FieldKeys]=DDpatternsAndKeys;
    %%
    dispFileStatus(DD.path)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DD=evalUserInput
   readfiles;
    setOutWindowIfNotUserSet;
    %----------------------------------------------------------------------
    function setOutWindowIfNotUserSet
        wesn={'west' 'east' 'south' 'north'};
        outmapchosen=isfield(DD.map.out,wesn);
        if ~all(outmapchosen)
            for sd=find(~outmapchosen);
                skydir=wesn{sd};
                DD.map.out.(skydir)=DD.map.in.(skydir);
            end
        end
    end
    %----------------------------------------------------------------------
    function readfiles
        B=INPUT;
        A=eval(['INPUT' B.template]);
        DD=mergeStruct2(A,B);
    end     
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pattern,FieldKeys]=DDpatternsAndKeys
    %% dir and file formats
    pattern.fname='CUT_yyyymmdd_SSSSsNNNNnWWWWwEEEEe.mat';
    pattern.prefix.cuts='CUT';
    pattern.prefix.conts='CONT';
    pattern.prefix.eddies='EDDIE';
    pattern.prefix.tracks='TRACK';
    %% fields that must end with .mean and .std - for output plot maps %
    FieldKeys.MeanStdFields= { ...
        'age';
        'dist.traj.fromBirth';
        'dist.traj.tillDeath';
        'dist.zonal.fromBirth';
        'dist.zonal.tillDeath';
        'dist.merid.fromBirth';
        'dist.merid.tillDeath';
        'radius.mean';
        'radius.zonal';
        'radius.meridional';
        'vel.traj';
        'vel.zonal';
        'vel.merid';
        'amp.to_contour';
        'amp.to_ellipse';
        'amp.to_mean';
        };
    %% fields 4 colorcoded track plots
    FieldKeys.trackPlots= { ...
        'isoper';
        'radius.mean';
        'radius.meridional';
        'radius.zonal';
        'age';
        'peak.amp.to_contour';
        'peak.amp.to_mean';
        'peak.amp.to_ellipse';
        };
    %% TODO
    FieldKeys.senses= { ...
        'AntiCycs';
        'Cycs';
        };
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=timestuff(T)
    T.from.num=datenum(T.from.str,'yyyymmdd');
    T.till.num=datenum(T.till.str,'yyyymmdd');
    T.span=T.till.num-T.from.num+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mkDirs(path)
    %%
    mkdirp(path.root);
    mkdirp(path.plots);
    mkdirp(path.code);
    mkdirp(path.codesubs);
    mkdirp(path.cuts.name);
    mkdirp(path.conts.name);
    mkdirp(path.eddies.name);
    mkdirp(path.tracks.name);
    mkdirp(path.analyzed.name);
    mkdirp(path.analyzedTracks.AC.name);
    mkdirp(path.analyzedTracks.C.name);
    mkdirp(path.Rossby.name);
    %%
    system(['cp ./*.m ' path.code]);
    system(['cp ./SUBS/*.m ' path.codesubs]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path=findfiles(DD)
    path=DD.path;
    path.root=['../data' path.OutDirBaseName '/'];
    path.plots=['../PLOTS/' path.OutDirBaseName '/'];
    path.code=[path.root, 'code/'];
    path.codesubs=[path.root, 'code/SUBS/'];
    path.cuts.name=[path.root, 'CUTS/'];
    path.conts.name=[path.root, 'CONTS/'];
    path.eddies.name=[path.root,'EDDIES/'];
    path.tracks.name=[path.root,'TRACKS/'];
    path.analyzed.name=[path.root,'ANALYZED/'];
    path.analyzedTracks.AC.name=[path.analyzed.name,'AntiCyclones/'];
    path.analyzedTracks.C.name=[path.analyzed.name,'Cyclones/'];
    path.Rossby.name=[path.root,'Rossby/'];
    path.Rossby.Nfile=[path.Rossby.name,'N.cdf'];
    %%
    mkDirs(path)
    %%
    [~,~,ext.raw]=fileparts(DD.map.in.fname); 
    patt=strsplit(DD.map.in.fname,'yyyymmdd');
    path.raw.files=dir([path.raw.name,patt{1},'*',patt{2}]);
    path.protoMaps.file=[path.root, 'protoMaps.mat'];
    path.meanU.file=[path.root, 'meanU.mat'];
    path.UV.files=dir([path.UV.name,'*.nc']);   
    path.cuts.files=dir([path.cuts.name,'*.mat']);
    path.conts.files=dir([path.conts.name,'*.mat']);
    path.eddies.files=dir([path.eddies.name,'*.mat']);
    path.tracks.files=dir([path.tracks.name,'*.mat']);
    path.analyzed.files=dir([path.analyzed.name,'*.mat']);
    path.analyzedTracks.AC.files=dir([path.analyzedTracks.AC.name,'*.mat']);
    path.analyzedTracks.C.files=dir([path.analyzedTracks.C.name,'*.mat']);
    path.Rossby.files=[dir([path.Rossby.name,'*.nc']); dir([path.Rossby.name,'*.mat'])];
    %%
    path.TempSalt.files=tempsalt(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function files=tempsalt(DD)
    try
        files=dir([DD.path.TempSalt.name,'*.nc']);
    catch  %#ok<CTCH>
        disp('found no Salt/Temp data. Skipping Rossby calcs..')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispFileStatus(p)
    disp(['found ' num2str(numel(p.raw.files)) ' files in ' p.raw.name])
    disp(['found ' num2str(numel(p.cuts.files)) ' files in ' p.cuts.name])
    disp(['found ' num2str(numel(p.conts.files)) ' files in ' p.conts.name])
    disp(['found ' num2str(numel(p.eddies.files)) ' files in ' p.eddies.name])
    disp(['found ' num2str(numel(p.tracks.files)) ' files in ' p.tracks.name])
    disp(['found ' num2str(numel(p.Rossby.files)) ' files in ' p.Rossby.name])
    disp(['found ' num2str(numel(p.analyzedTracks.AC.files)) ' files in ' p.analyzedTracks.AC.name])
    disp(['found ' num2str(numel(p.analyzedTracks.C.files)) ' files in ' p.analyzedTracks.C.name])
    disp(['found ' num2str(numel(p.analyzed.files)) ' files in ' p.analyzed.name])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
