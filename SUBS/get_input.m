function [DD]=get_input
	%%
	DD=input_vars;
	%%
	DD.time.from.num=datenum(DD.time.from.str,'yyyymmdd');
	DD.time.till.num=datenum(DD.time.till.str,'yyyymmdd');
	DD.time.span=DD.time.till.num-DD.time.from.num+1;
	%%
<<<<<<< HEAD
	DD.path.cuts.name=[DD.path.root, 'CUTS/'];
=======
	DD.path.code=[DD.path.root, 'code/'];
    DD.path.codesubs=[DD.path.root, 'code/SUBS/'];
    DD.path.cuts.name=[DD.path.root, 'CUTS/'];
>>>>>>> master
	DD.path.conts.name=[DD.path.root, 'CONTS/'];
	DD.path.eddies.name=[DD.path.root,'EDDIES/'];
	DD.path.tracks.name=[DD.path.root,'TRACKS/'];
	DD.path.analyzed.name=[DD.path.root,'ANALYZED/'];
	DD.path.analyzedTracks.AC.name=[DD.path.analyzed.name,'AntiCyclones/'];
	DD.path.analyzedTracks.C.name=[DD.path.analyzed.name,'Cyclones/'];
    DD.path.Rossby.name=[DD.path.root,'Rossby/'];
	
   	%%	
    mkdirp(DD.path.code);
    mkdirp(DD.path.codesubs);
	mkdirp(DD.path.cuts.name);
	mkdirp(DD.path.conts.name);
	mkdirp(DD.path.eddies.name);
	mkdirp(DD.path.tracks.name);
	mkdirp(DD.path.analyzed.name);	
	mkdirp(DD.path.analyzedTracks.AC.name);	
	mkdirp(DD.path.analyzedTracks.C.name);	
    mkdirp(DD.path.Rossby.name);	
    %%
    system(['cp ./*.m ' DD.path.code]);
    system(['cp ./SUBS/*.m ' DD.path.codesubs]);
	%%
	 DD.path.protoMaps.file=[DD.path.root, 'protoMaps.mat'];
	 DD.path.meanU.file=[DD.path.root, 'meanU.mat'];
	DD.path.raw.files=dir([DD.path.raw.name,'*.nc']);
	DD.path.cuts.files=dir([DD.path.cuts.name,'*.mat']);
	DD.path.conts.files=dir([DD.path.conts.name,'*.mat']);
	DD.path.eddies.files=dir([DD.path.eddies.name,'*.mat']);
	DD.path.tracks.files=dir([DD.path.tracks.name,'*.mat']);
	DD.path.analyzed.files=dir([DD.path.analyzed.name,'*.mat']);
	DD.path.analyzedTracks.AC.files=dir([DD.path.analyzedTracks.AC.name,'*.mat']);
	DD.path.analyzedTracks.C.files=dir([DD.path.analyzedTracks.C.name,'*.mat']);
    DD.path.Rossby.files=dir([DD.path.Rossby.name,'*.nc']);
	%%
	DD.pattern.in='CUT_yyyymmdd_SSSSsNNNNnWWWWwEEEEe.mat'; 
	try
		DD.path.TempSalt.files=dir([DD.path.TempSalt.name,'*.nc']);
	catch  %#ok<CTCH>
		disp('found no Salt/Temp data. Skipping Rossby calcs..')
	end
end
