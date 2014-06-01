function [DD]=get_input
    %%
    DD=input_vars;
    %%
    DD.time=catstruct(DD.time, timestuff(DD.time));
    %%
    
    DD.path=catstruct(DD.path,findfiles(DD.path));
    %%
    DD.pattern.in='CUT_yyyymmdd_SSSSsNNNNnWWWWwEEEEe.mat';  
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
function mkDirs(Path)
    %%
    mkdirp(Path.code);
    mkdirp(Path.codesubs);
    mkdirp(Path.cuts.name);
    mkdirp(Path.conts.name);
    mkdirp(Path.eddies.name);
    mkdirp(Path.tracks.name);
    mkdirp(Path.analyzed.name);
    mkdirp(Path.analyzedTracks.AC.name);
    mkdirp(Path.analyzedTracks.C.name);
    mkdirp(Path.Rossby.name);
    %%
    system(['cp ./*.m ' Path.code]);
    system(['cp ./SUBS/*.m ' Path.codesubs]);
end

function Path=findfiles(Path)
    Path.code=[Path.root, 'code/'];
    Path.codesubs=[Path.root, 'code/SUBS/'];
    Path.cuts.name=[Path.root, 'CUTS/'];
    Path.conts.name=[Path.root, 'CONTS/'];
    Path.eddies.name=[Path.root,'EDDIES/'];
    Path.tracks.name=[Path.root,'TRACKS/'];
    Path.analyzed.name=[Path.root,'ANALYZED/'];
    Path.analyzedTracks.AC.name=[Path.analyzed.name,'AntiCyclones/'];
    Path.analyzedTracks.C.name=[Path.analyzed.name,'Cyclones/'];
    Path.Rossby.name=[Path.root,'Rossby/'];
    %%
    mkDirs(Path)
    %%
    Path.protoMaps.file=[Path.root, 'protoMaps.mat'];
    Path.meanU.file=[Path.root, 'meanU.mat'];
    Path.raw.files=dir([Path.raw.name,'*.nc']);
    Path.cuts.files=dir([Path.cuts.name,'*.mat']);
    Path.conts.files=dir([Path.conts.name,'*.mat']);
    Path.eddies.files=dir([Path.eddies.name,'*.mat']);
    Path.tracks.files=dir([Path.tracks.name,'*.mat']);
    Path.analyzed.files=dir([Path.analyzed.name,'*.mat']);
    Path.analyzedTracks.AC.files=dir([Path.analyzedTracks.AC.name,'*.mat']);
    Path.analyzedTracks.C.files=dir([Path.analyzedTracks.C.name,'*.mat']);
    Path.Rossby.files=dir([Path.Rossby.name,'*.nc']);
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
