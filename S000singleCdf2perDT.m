%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 25-Apr-2014 12:01:45
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Madeleine Version
function S000singleCdf2perDT
	%% init dependencies
	addpath(genpath('./'));
	%% get user input
	DD = initialise;
	%% get madeleine's data
	[raw]=cdfData(DD);
	%% get geo stuff
	[DD,raw]=geostuff(raw,DD);
	%% thread distro
	disp('working through all timesteps for now!')
	DD.threads.lims=thread_distro(DD.threads.num,numel(raw.TIME));
	%% start threads
	init_threads(DD.threads.num);
	%% spmd
	main(DD,raw)
	%% save brunt väisälä
	saveN(DD,raw)
	%% save info
	save_info(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveN(DD,raw)
	N=double(squeeze(nc_varget(raw.file.in,DD.map.in.keys.N,[0 0 0 0],[1 inf inf inf])));
	Nfile=[DD.path.Rossby.name 'N.cdf'];
	NCoverwriteornot(Nfile);
	nc_adddim(Nfile,'i_index',DD.map.window.size.X);
	nc_adddim(Nfile,'j_index',DD.map.window.size.Y);
	nc_adddim(Nfile,'k_index',DD.map.window.size.Z);
	%% N
	varstruct.Name = DD.map.in.keys.N;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'k_index','j_index','i_index' };
	nc_addvar(Nfile,varstruct);
	%% put
	nc_varput(Nfile,varstruct.Name,N);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,raw)
	if DD.debugmode
		spmd_body(DD,raw);
	else
		spmd(DD.threads.num)
			spmd_body(DD,raw);
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD,raw)
	CC=(DD.threads.lims(labindex,1):DD.threads.lims(labindex,2));
	%% loop over files
	[T]=disp_progress('init','preparing raw data');
	for cc=CC
		[T]=disp_progress('calc',T,numel(CC),100);
		%% get current SSH
		raw.grids.ssh=squeeze(nc_varget(raw.file.in,DD.map.in.keys.ssh,[cc-1,raw.SSHzIdx-1,0,0],[1,1,inf,inf]));
		operateDay(raw,DD,cc);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [raw]=cdfData(DD)
	raw.file.in=[DD.path.raw.name	,DD.map.in.cdfName];
% 	raw.info=ncInfoAll(raw.file.in);
% 	for info=fieldnames(raw.info)'; disp(raw.info.(info{1})); end
	disp(['setting user start date - ' DD.time.from.str ' - as start date!'])
	startTime=DD.time.from.num;
	keys=DD.map.in.keys;
	raw.TIME=nc_varget(raw.file.in,keys.time);
	raw.TIME=raw.TIME-raw.TIME(1)+startTime;
	raw.XT=nc_varget(raw.file.in,keys.x);
	raw.YT=nc_varget(raw.file.in,keys.y);
	raw.ZT=nc_varget(raw.file.in,keys.z);
	[~,raw.SSHzIdx]=min(abs(raw.ZT-DD.parameters.SSHAdepth));
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,raw]=geostuff(raw,DD)
	[raw.grids.XX,raw.grids.YY]=meshgrid(raw.XT,raw.YT);
	raw.grids.lat=rad2deg(raw.grids.YY./earthRadius) + DD.parameters.boxlims.south;
	raw.grids.lon=rad2deg(raw.grids.XX./(cosd(raw.grids.lat)*earthRadius)) +  DD.parameters.boxlims.west;
	if max(diff(raw.grids.lon(:)))>300, error('dont put window on -180/180 meridian!'); end
	[raw.grids.DY,raw.grids.DX]=DYDX(raw.grids.lat,raw.grids.lon);
	%% reset to exact values
	DD.map.in.west=min(raw.grids.lon(:));
	DD.map.in.east=max(raw.grids.lon(:));
	DD.map.in.south=min(raw.grids.lat(:));
	DD.map.in.north=max(raw.grids.lat(:));
	%% reset out maps
	DD.map.out=getOutMapRes(DD.map.out);
	DD.map.out.west=DD.map.in.west;
	DD.map.out.east=DD.map.in.east;
	DD.map.out.south=DD.map.in.south;
	DD.map.out.north=DD.map.in.north;
	%% use full map
	[Y,X]=size(raw.grids.lon);
	DD.map.window.size.X=X;
	DD.map.window.size.Y=Y;
	DD.map.window.limits.west=1;
	DD.map.window.limits.east=X;
	DD.map.window.limits.south=1;
	DD.map.window.limits.north=Y;
	DD.map.window.size.Z=numel(raw.ZT);
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
function saveraw(DD,raw)
	NCoverwriteornot(raw.file.out);
	nc_adddim(raw.file.out,'i_index',DD.map.window.size.X);
	nc_adddim(raw.file.out,'j_index',DD.map.window.size.Y);
	%% lat
	varstruct.Name = DD.map.in.keys.lat;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'j_index','i_index' };
	nc_addvar(raw.file.out,varstruct);
	%% lon
	varstruct.Name = DD.map.in.keys.lon;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'j_index','i_index' };
	nc_addvar(raw.file.out,varstruct);
	%% ssh
	varstruct.Name = DD.map.in.keys.ssh;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'j_index','i_index' };
	nc_addvar(raw.file.out,varstruct);
	%%----------put-----------------
	%%------------------------------
	nc_varput(raw.file.out,DD.map.in.keys.lat,raw.grids.lat);
	nc_varput(raw.file.out,DD.map.in.keys.lon,raw.grids.lon);
	nc_varput(raw.file.out,DD.map.in.keys.ssh,raw.grids.ssh);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function operateDay(raw,DD,cc)
	%% set up output file
	tt=raw.TIME(cc);
	timestr=datestr(tt,'yyyymmdd');
	path=DD.path.raw.name;
	fo=DD.map.in.fname;
	fo=strrep(fo,'yyyymmdd',timestr);
	raw.file.out=[path, fo];
	if exist(raw.file.out,'file'), return; end
	%%
	foulIdx=(raw.grids.ssh>1000 | raw.grids.ssh<-1000 | isnan(raw.grids.ssh));
	raw.grids.ssh=double(NeighbourValue(foulIdx, raw.grids.ssh));
	%%
	saveraw(DD,raw);
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
