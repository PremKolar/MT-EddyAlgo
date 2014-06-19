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
    [RAW]=cdfData(DD);
    %% get geo stuff
    [DD,RAW]=geostuff(RAW,DD);
    %% thread distro
    disp('working through all timesteps for now!')
    DD.threads.lims=thread_distro(DD.threads.num,numel(RAW.TIME));
    %% start threads
    init_threads(DD.threads.num);
    %% spmd
    main(DD,RAW)
    %% save info
    save_info(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,RAW)
    if DD.debugmode
        spmd_body(DD,RAW);
    else
        spmd(DD.threads.num)
            spmd_body(DD,RAW);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD,RAW)
    CC=(DD.threads.lims(labindex,1):DD.threads.lims(labindex,2));
    %% loop over files
    [T]=disp_progress('init','preparing raw data');
    for cc=CC
        [T]=disp_progress('calc',T,numel(CC),100);
        %% get current SSH 
        RAW.grids.ssh=squeeze(nc_varget(RAW.file.in,'SSHA',[cc-1,RAW.SSHzIdx-1,0,0],[1,1,inf,inf]));
        operateDay(RAW,DD,cc);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RAW]=cdfData(DD)
    RAW.file.in=[DD.path.raw.name	,DD.map.in.cdfName	];
    RAW.info=ncInfoAll(RAW.file.in);
    for info=fieldnames(RAW.info)'; disp(RAW.info.(info{1})); end
    disp(['setting user start date - ' DD.time.from.str ' - as start date!'])
    startTime=DD.time.from.num;
    RAW.TIME=nc_varget(RAW.file.in,'TIME');
    RAW.TIME=RAW.TIME-RAW.TIME(1)+startTime;
    RAW.XT=nc_varget(RAW.file.in,'XT');
    RAW.YT=nc_varget(RAW.file.in,'YT');
    RAW.ZT=nc_varget(RAW.file.in,'ZT');    
   [~,RAW.SSHzIdx]=min(abs(RAW.ZT-DD.parameters.SSHAdepth));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,RAW]=geostuff(RAW,DD)
    [RAW.grids.XX,RAW.grids.YY]=meshgrid(RAW.XT,RAW.YT);
    RAW.grids.lat=rad2deg(RAW.grids.YY./earthRadius) + DD.parameters.boxlims.south;
    RAW.grids.lon=rad2deg(RAW.grids.XX./(cosd(RAW.grids.lat)*earthRadius)) +  DD.parameters.boxlims.west;
    if max(diff(RAW.grids.lon(:)))>300, error('dont put window on -180/180 meridian!'); end
    [RAW.grids.DY,RAW.grids.DX]=DYDX(RAW.grids.lat,RAW.grids.lon);
    %% reset to exact values
    DD.map.in.west=min(RAW.grids.lon(:));
    DD.map.in.east=max(RAW.grids.lon(:));
    DD.map.in.south=min(RAW.grids.lat(:));
    DD.map.in.north=max(RAW.grids.lat(:));
    %% reset out maps
    DD.map.out=getOutMapRes(DD.map.out);
    DD.map.out.west=DD.map.in.west;
    DD.map.out.east=DD.map.in.east;
    DD.map.out.south=DD.map.in.south;
    DD.map.out.north=DD.map.in.north;
    %% use full map
    [Y,X]=size(RAW.grids.lon);
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
function saveRAW(DD,RAW)
    NCoverwriteornot(RAW.file.out);
    nc_adddim(RAW.file.out,'i_index',DD.map.window.size.X);
    nc_adddim(RAW.file.out,'j_index',DD.map.window.size.Y);
    %% lat
    varstruct.Name = DD.map.in.keys.lat;
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'j_index','i_index' };
    nc_addvar(RAW.file.out,varstruct);
    %% lon
    varstruct.Name = DD.map.in.keys.lon;
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'j_index','i_index' };
    nc_addvar(RAW.file.out,varstruct);
    %% ssh
    varstruct.Name = DD.map.in.keys.ssh;
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'j_index','i_index' };
    nc_addvar(RAW.file.out,varstruct);
    %%----------put-----------------
    %%------------------------------
    nc_varput(RAW.file.out,DD.map.in.keys.lat,RAW.grids.lat);
    nc_varput(RAW.file.out,DD.map.in.keys.lon,RAW.grids.lon);
    nc_varput(RAW.file.out,DD.map.in.keys.ssh,RAW.grids.ssh);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function operateDay(RAW,DD,cc)
    %% set up output file
    tt=RAW.TIME(cc);
    timestr=datestr(tt,'yyyymmdd');
    path=DD.path.raw.name;
    fo='RAWyyyymmdd.nc';
    fo=strrep(fo,'yyyymmdd',timestr);
    RAW.file.out=[path, fo];
    if exist(RAW.file.out,'file'), return; end
    %%
    RAW.grids.ssh(RAW.grids.ssh>1000)=nan;RAW.grids.ssh(RAW.grids.ssh<-1000)=nan;
    RAW.grids.ssh=double(RAW.grids.ssh/DD.map.in.ssh_unitFactor);
    %%
    saveRAW(DD,RAW);
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
