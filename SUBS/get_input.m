function [DD]=get_input
    %%
    DD=input_vars;
    %%
    DD.time=catstruct(DD.time, timestuff(DD.time));
    %%
    DD.path=catstruct(DD.path,findfiles(DD));
    %%
    DD.pattern.fname='CUT_yyyymmdd_SSSSsNNNNnWWWWwEEEEe.mat';
    DD.pattern.prefix.cuts='CUT';
    DD.pattern.prefix.conts='CONT';
    DD.pattern.prefix.eddies='EDDIE';
    DD.pattern.prefix.tracks='TRACK';
    %%
    DD.path.TempSalt.files=tempsalt(DD);
    %%
    dispFileStatus(DD.path)    
end
function T=timestuff(T)
    T.from.num=datenum(T.from.str,'yyyymmdd');
    T.till.num=datenum(T.till.str,'yyyymmdd');
    T.span=T.till.num-T.from.num+1;
end
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
    %%
    mkDirs(path)
    %%    
    [~,~,ext.raw]=fileparts(DD.map.in.fname);    
    path.protoMaps.file=[path.root, 'protoMaps.mat'];
    path.meanU.file=[path.root, 'meanU.mat'];
    path.raw.files=dir([path.raw.name,'*',ext.raw]);
    path.cuts.files=dir([path.cuts.name,'*.mat']);
    path.conts.files=dir([path.conts.name,'*.mat']);
    path.eddies.files=dir([path.eddies.name,'*.mat']);
    path.tracks.files=dir([path.tracks.name,'*.mat']);
    path.analyzed.files=dir([path.analyzed.name,'*.mat']);
    path.analyzedTracks.AC.files=dir([path.analyzedTracks.AC.name,'*.mat']);
    path.analyzedTracks.C.files=dir([path.analyzedTracks.C.name,'*.mat']);
    path.Rossby.files=dir([path.Rossby.name,'*.nc']);
end
function files=tempsalt(DD)
    try
        files=dir([DD.path.TempSalt.name,'*.nc']);
    catch  %#ok<CTCH>
        disp('found no Salt/Temp data. Skipping Rossby calcs..')
    end
end
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
