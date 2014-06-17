%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 25-Apr-2014 12:01:45
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Madeleine Version
function S00_prep_data
	%% init dependencies
	addpath(genpath('./'));
	%% get user input
	DD = initialise;
	%% get madeleine's data
	[MadFile,MF]=madsData(DD);
	%% get geo stuff
	[DD,MF]=geostuff(MF,DD);
	%% thread distro
	DD.threads.lims=thread_distro(DD.threads.num,numel(MF.TIME));
	%% start threads
	init_threads(DD.threads.num);
	%% create out dir
	[~,~]=mkdir(DD.path.cuts.name);
	%% spmd
	spmd
	spmdBlock(DD,MadFile,MF);
	end
	%% save info
	save_info(DD);
end
function spmdBlock(DD,MadFile,MF)
		CC=(DD.threads.lims(labindex,1):DD.threads.lims(labindex,2));
		E325=nc_varget(MadFile,'E325');
		%% loop over files
		[T]=disp_progress('init','preparing raw data');
		for cc=CC			
			[T]=disp_progress('calc',T,numel(CC),4242);
			operateDay(squeeze(E325(cc,:,:)),MF,DD,cc);
		end	
end
function [MadFile,MF]=madsData(DD)
	MadFile=[DD.path.root	,'psvar.cdf'];
	nc_info(MadFile);
	MF.TIME=nc_varget(MadFile,'TIME')+datenum('1900','yyyy');
	MF.XT_bnds=nc_varget(MadFile,'XT_bnds');
	MF.YT_bnds=nc_varget(MadFile,'YT_bnds');
	MF.XT=nc_varget(MadFile,'XT');
	MF.YT=nc_varget(MadFile,'YT');
end
function [DD,MF]=geostuff(MF,DD)
	[MF.grids.XX,MF.grids.YY]=meshgrid(MF.XT,MF.YT);
	MF.grids.lat=rad2deg(MF.grids.YY./earthRadius);
	MF.grids.lon=rad2deg(MF.grids.XX./(cosd(MF.grids.lat)*earthRadius));
	[MF.grids.DY,MF.grids.DX]=DYDX(MF.grids.lat,MF.grids.lon);
	DD.map.west=min(MF.grids.lon(:));
	DD.map.east=max(MF.grids.lon(:));
	DD.map.south=min(MF.grids.lat(:));
	DD.map.north=max(MF.grids.lat(:));
	[Y,X]=size(MF.grids.lon);
	DD.map.window.size.X=X;
	DD.map.window.size.Y=Y;
	DD.map.window.limits.west=1;
	DD.map.window.limits.east=X;
	DD.map.window.limits.south=1;
	DD.map.window.limits.north=Y;
end
function operateDay(SSH,MF,DD,cc)
	%% set up output file
	tt=MF.TIME(cc);
	timestr=datestr(tt,'yyyymmdd');
	path=DD.path.cuts.name;
	geo=DD.map;
	file.out=strrep(DD.pattern.in	,'SSSS',sprintf('%04d',round(geo.south)) );
	file.out=strrep(file.out, 'NNNN',sprintf('%04d',round(geo.north) ));
	file.out=strrep(file.out, 'WWWW',sprintf('%04d',round(geo.west) ));
	file.out=strrep(file.out, 'EEEE',sprintf('%04d',round(geo.east)) );
	file.out=[path, strrep(file.out, 'yyyymmdd',timestr)];
	 if exist(file.out,'file'), return; end
	%%
	SSH(SSH>10000)=nan;
	SSH(SSH<-10000)=nan;
	MF.grids.ssh=double(SSH/DD.map.ssh_unitFactor);
	%%
	MF.params.full_globe.x=false;
	%%
	save(file.out,'-struct','MF')
end
function [DY,DX]=DYDX(LAT,LON)
	%% grid increment sizes
	DY=deg2rad(abs(diff(double(LAT),1,1)))*earthRadius;
	DX=deg2rad(abs(diff(double(LON),1,2)))*earthRadius.*cosd(LAT(:,1:end-1));
	%% append one line/row to have identical size as other fields
	DY=DY([1:end,end],:);
	DX=DX(:,[1:end,end]);
	%% correct 360° crossings
	seamcrossflag=DX>100*median(DX(:));
	DX(seamcrossflag)=abs(DX(seamcrossflag) - 2*pi*earthRadius.*cosd(LAT(seamcrossflag)));
end
