%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 25-Apr-2014 12:01:45
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Madeleine Version
function S000cdf2daily
	%% init dependencies
	addpath(genpath('./'));
	%% get user input
	DD = initialise;
	%% get madeleine's data
	[MF]=cdfData(DD);
	%% get geo stuff
	[DD,MF]=geostuff(MF,DD);
	%% thread distro
	disp('working through all timesteps for now!')
	DD.threads.lims=thread_distro(DD.threads.num,numel(MF.TIME));
	%% start threads
	init_threads(DD.threads.num);
	%% spmd
	main(DD,MF)
	%% save info
	save_info(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,MF)
	if DD.debugmode
		spmd_body(DD,MF);
	else
		spmd(DD.threads.num)
			spmd_body(DD,MF);
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD,MF)
	CC=(DD.threads.lims(labindex,1):DD.threads.lims(labindex,2));
	%% loop over files
	[T]=disp_progress('init','preparing raw data');
	for cc=CC
		MF.E325=squeeze(nc_varget(MF.file.in,'E325',[cc-1,0,0],[1,inf,inf]));
		[T]=disp_progress('calc',T,numel(CC),100);
		operateDay(MF,DD,cc);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MF]=cdfData(DD)
	MF.file.in=[DD.path.raw.name	,DD.map.in.cdfName	];
	MF.nc_info=nc_info(MF.file.in);
	disp(['setting user start date - ' DD.time.from.str ' - as start date!'])
	startTime=DD.time.from.num;
	MF.TIME=nc_varget(MF.file.in,'TIME');
	MF.TIME=MF.TIME-MF.TIME(1)+startTime;
	MF.XT_bnds=nc_varget(MF.file.in,'XT_bnds');
	MF.YT_bnds=nc_varget(MF.file.in,'YT_bnds');
	MF.XT=nc_varget(MF.file.in,'XT');
	MF.YT=nc_varget(MF.file.in,'YT');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,MF]=geostuff(MF,DD)
	[MF.grids.XX,MF.grids.YY]=meshgrid(MF.XT,MF.YT);
	MF.grids.lat=rad2deg(MF.grids.YY./earthRadius) + DD.parameters.boxlims.south;
	MF.grids.lon=rad2deg(MF.grids.XX./(cosd(MF.grids.lat)*earthRadius)) +  DD.parameters.boxlims.west;
	if max(diff(MF.grids.lon(:)))>300, error('dont put window on -180/180 meridian!'); end
	[MF.grids.DY,MF.grids.DX]=DYDX(MF.grids.lat,MF.grids.lon);
	%% reset to exact values
	DD.map.in.west=min(MF.grids.lon(:));
	DD.map.in.east=max(MF.grids.lon(:));
	DD.map.in.south=min(MF.grids.lat(:));
	DD.map.in.north=max(MF.grids.lat(:));
	%% reset out maps
	DD.map.out=getOutMapRes(DD.map.out);
	DD.map.out.west=DD.map.in.west;
	DD.map.out.east=DD.map.in.east;
	DD.map.out.south=DD.map.in.south;
	DD.map.out.north=DD.map.in.north;
	%% use full map
	[Y,X]=size(MF.grids.lon);
	DD.map.window.size.X=X;
	DD.map.window.size.Y=Y;
	DD.map.window.limits.west=1;
	DD.map.window.limits.east=X;
	DD.map.window.limits.south=1;
	DD.map.window.limits.north=Y;
	%% info
	mapInfo(Y,X,DD.map.in)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=getOutMapRes(map)
	londiff=(map.east - map.west );
	latdiff=(map.east - map.west );
	res.x=(map.X-1)/londiff;
	res.y=(map.Y-1)/latdiff;
	N.X=londiff*res.x +1;
	N.Y=latdiff*res.y +1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mapInfo(Y,X,map)
	sprintf('built %ix%i grid',Y,X);
	sprintf('spanning %05.1fW:%05.1fE / %05.1fS:%05.1fN',map.west,map.east,map.south,map.north);
	sleep(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveMF(DD,MF)
	NCoverwriteornot(MF.file.out);
	nc_adddim(MF.file.out,'i_index',DD.map.window.size.X);
	nc_adddim(MF.file.out,'j_index',DD.map.window.size.Y);
	%% lat
	varstruct.Name = DD.map.in.keys.lat;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'j_index','i_index' };
	nc_addvar(MF.file.out,varstruct);
	%% lon
	varstruct.Name = DD.map.in.keys.lon;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'j_index','i_index' };
	nc_addvar(MF.file.out,varstruct);
	%% ssh
	varstruct.Name = DD.map.in.keys.ssh;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'j_index','i_index' };
	nc_addvar(MF.file.out,varstruct);
	%%----------put-----------------
	%%------------------------------
	nc_varput(MF.file.out,DD.map.in.keys.lat,MF.grids.lat);
	nc_varput(MF.file.out,DD.map.in.keys.lon,MF.grids.lon);
	nc_varput(MF.file.out,DD.map.in.keys.ssh,MF.grids.ssh);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function operateDay(MF,DD,cc)
	%% set up output file
	tt=MF.TIME(cc);
	timestr=datestr(tt,'yyyymmdd');
	path=DD.path.raw.name;
	fo='RAWyyyymmdd.nc';
	fo=strrep(fo,'yyyymmdd',timestr);
	MF.file.out=[path, fo];
	if exist(MF.file.out,'file'), return; end
	%%
	MF.E325(MF.E325>10000)=nan;
	MF.E325(MF.E325<-10000)=nan;
	MF.grids.ssh=double(MF.E325/DD.map.in.ssh_unitFactor);
	%%
	MF.params.full_globe.x=false;
	%%
	saveMF(DD,MF);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
