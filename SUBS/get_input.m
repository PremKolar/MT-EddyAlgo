function [DD]=get_input
    sprintf(['\n getting user input...']);
    %%
    DD=evalUserInput;
    %%
    DD.time=catstruct(DD.time, timestuff(DD.time));
    %%
    sprintf(['\n setting internal parameters...']);
    [DD.pattern, DD.FieldKeys]=DDpatternsAndKeys;
    %%
    sprintf(['\n scanning data...']);
    DD.path=catstruct(DD.path,findfiles(DD));
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
function mkDirs(path,senses)
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
    mkdirp(path.analyzedTracks.(senses{1}).name);
    mkdirp(path.analyzedTracks.(senses{2}).name);
    mkdirp(path.Rossby.name);
    %%
    system(['cp ./*.m ' path.code]);
    system(['cp ./SUBS/*.m ' path.codesubs]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PATH=findfiles(DD)
    senses=DD.FieldKeys.senses;
    
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
    PATH.analyzedTracks.(senses{1}).name=[PATH.analyzed.name,senses{1},'/'];
    PATH.analyzedTracks.(senses{2}).name=[PATH.analyzed.name,senses{2},'/'];
    PATH.analyzedTracks.C.name=[PATH.analyzed.name,'Cyclones/'];
    PATH.Rossby.name=[PATH.root,'Rossby/'];
    PATH.Rossby.Nfile=[PATH.Rossby.name,'N.cdf'];
    %%
    mkDirs(PATH,senses)
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
    PATH.analyzedTracks.(senses{1}).files=dir([PATH.analyzedTracks.(senses{1}).name,'*.mat']);
    PATH.analyzedTracks.(senses{2}).files=dir([PATH.analyzedTracks.(senses{2}).name,'*.mat']);
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
    FN=fieldnames(p)';
    for ii=1:numel(FN);fn=FN{ii};
        if isfield(p.(fn),'files') && isfield(p.(fn),'name')
            disp(['found ' num2str(numel(p.(fn).files)) ' files in ' p.(fn).name]);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
