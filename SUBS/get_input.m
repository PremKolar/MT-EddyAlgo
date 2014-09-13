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
    %% Rossby
    FieldKeys.Rossby = { ...
    'RossbyPhaseSpeed'   ;
    'RossbyRadius' ;
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
function PATH=findfiles(DD)
    PATH=DD.path;
    PATH.root=['../data' PATH.OutDirBaseName '/'];
    PATH.plots=['../PLOTS/' PATH.OutDirBaseName '/'];
    PATH.code=[PATH.root, 'code/'];
    PATH.codesubs=[PATH.root, 'code/SUBS/'];
    PATH.cuts.name=[PATH.root, 'CUTS/'];
    PATH.conts.name=[PATH.root, 'CONTS/'];
    PATH.eddies.name=[PATH.root,'EDDIES/'];
    PATH.tracks.name=[PATH.root,'TRACKS/'];
    PATH.analyzed.name=[PATH.root,'ANALYZED/'];
    PATH.analyzedTracks.AC.name=[PATH.analyzed.name,'AntiCyclones/'];
    PATH.analyzedTracks.C.name=[PATH.analyzed.name,'Cyclones/'];
    PATH.Rossby.name=[PATH.root,'Rossby/'];
    PATH.Rossby.Nfile=[PATH.Rossby.name,'N.cdf'];
    %%
    mkDirs(PATH)
    %%
    [~,~,ext.raw]=fileparts(DD.map.in.fname);
    patt=strsplit(DD.map.in.fname,'yyyymmdd');
    PATH.raw.files=dir([PATH.raw.name,patt{1},'*',ext.raw]);
    PATH.protoMaps.file=[PATH.root, 'protoMaps.mat'];
    PATH.meanU.file=[PATH.root, 'meanU.mat'];
    PATH.UV.files=dir([PATH.UV.name,'*.nc']);
    PATH.cuts.files=dir([PATH.cuts.name,'*.mat']);
    PATH.conts.files=dir([PATH.conts.name,'*.mat']);
    PATH.eddies.files=dir([PATH.eddies.name,'*.mat']);
    PATH.tracks.files=dir([PATH.tracks.name,'*.mat']);
    PATH.analyzed.files=dir([PATH.analyzed.name,'*.mat']);
    PATH.analyzedTracks.AC.files=dir([PATH.analyzedTracks.AC.name,'*.mat']);
    PATH.analyzedTracks.C.files=dir([PATH.analyzedTracks.C.name,'*.mat']);
    PATH.Rossby.files=[dir([PATH.Rossby.name,'*.nc']); dir([PATH.Rossby.name,'*.mat'])];
    %%
    PATH.TempSalt.files=tempsalt(DD);
    
    PATH.windowFile=[PATH.root 'window.mat'];
  
    
    
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
