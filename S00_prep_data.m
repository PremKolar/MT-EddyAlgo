%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare SSH data
% reads user input from input_vars.m and map_vars.m
% This is the only file that would have to adapted to the strcture of the
% input SSH data. This step is not officially part of the program. Use a
% copy of this and adapt to your data, so that S01 gets the required input
% structure.
function S00_prep_data
    %% set up
    [DD]=set_up;
    %% spmd
    main(DD)
    %% save info
    save_info(DD)
end

function main(DD)
    if DD.debugmode
        spmd_body(DD);
    else
        spmd(DD.threads.num)
            spmd_body(DD);
        end
    end
end


function [DD]=set_up
    %% init dependencies
    addpath(genpath('./'));
    %% get user input
    DD = initialise('raw');
    %% get user map input
    DD.map=map_vars;
    %% get sample window
    [DD.map.window]=GetWindow(SampleFile(DD),DD);
end

function spmd_body(DD)
    %% distro chunks to threads
     [TT]=SetThreadVar(DD);   
    %% loop over files
    [T]=disp_progress('init','preparing raw data');
    for cc=1:numel(TT);
        %%
        [T]=disp_progress('calc',T,numel(TT),100);
        %% get data
        file=GetCurrentFile(TT(cc),DD)  ;     
        %% cut data
        [CUT,readable]=CutMap(file,DD); if ~readable, continue; end
        %% write data
        WriteFileOut(file.out,CUT);
    end
    
end
%% window functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window,readable]=GetWindow(file,DD)
    disp('assuming identical LON/LAT for all files!!!')
    %% get data
    [grids,readable]=GetFields(file.in,DD.map.pattern); if ~readable, return; end
    %% find window mask
    window=FindWindowMask(grids,DD.map.geo);
    %% find rectangle enclosing all applicable data
    window.limits=FindRectangle(window.flag);
    %% size
    window.size=WriteSize(window.limits);
end
function S=WriteSize(lims)
    S.X = lims.east-lims.west   +1;
    S.Y = lims.north-lims.south +1;
end
function [F,readable]=GetFields(file,keys)
    F=struct;
    readable=true;  
    try
        F.LON = CorrectLongitude(nc_varget(file,keys.lon)); %TODO: from user input AND: flexible varget function
        F.LAT = nc_varget(file,keys.lat);
        F.SSH = squeeze(nc_varget(file,keys.ssh));
        if numel(F.LON)~=numel(F.SSH)
          F.LON=repmat(F.LON',size(F.SSH,1),1);
           F.LAT=repmat(F.LAT,1,size(F.SSH,2));
        end
    catch void
        readable=false;
        warning(void.identifier,	['cant read ',file,', skipping!'])
        disp(void.message)
    end
end
function [LON]=CorrectLongitude(LON)
    % longitude(-180:180) concept is to be used!
    if max(LON(:))>180
        lontrans=true;
    else
        lontrans=false;
    end
    if lontrans
        LON(LON>180)=LON(LON>180)-360;
    end
end
function window=FindWindowMask(F,M)
    %% tag all grid points fullfilling all desired lat/lon limits
    if M.east>M.west
        window.flag= F.LON>=M.west & F.LON<=M.east & F.LAT>=M.south & F.LAT<=M.north ;
    elseif M.west>M.east  %crossing 180 meridian
        window.flag=((F.LON>=M.west & F.LON<=180) | (F.LON>=-180 & F.LON<=M.east)) & F.LAT>=M.south & F.LAT<=M.north ;
    end
end
function limits=FindRectangle(flag)
    %% find index limits
    [rows,cols]=find(flag);
    limits.west=min(cols);
    limits.east=max(cols);
    limits.north=max(rows);
    limits.south=min(rows);
end
%% Cutting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CUT,readable]=CutMap(file,DD)
    CUT=struct;
    %% get data
    [raw_fields,readable]=GetFields(file.in,DD.map.pattern); if ~readable, return; end
    %% cut
    [CUT]=SeamOrGlobe(raw_fields,DD.map.window);
    %% nan out land and make SI
    CUT.grids.SSH=nanLand(CUT.grids.SSH,DD.map.SSH_unitFactor);
    %% get distance fields
    [CUT.grids.DY,CUT.grids.DX]=DYDX(CUT.grids.LAT,CUT.grids.LON);
end
function [OUT]=SeamOrGlobe(IN,window)
    %% shorthands
    Wflag=window.flag;
    Wlin=window.limits;
    %% full globe?
    full_globe.x=false;
    full_globe.y=false;
    if (Wlin.east-Wlin.west+1==size(Wflag,2))
        full_globe.x=true;
        if (Wlin.north-Wlin.south+1==size(Wflag,1))
            full_globe.y=true;
        end
        OUT=AppendIfFullZonal(IN,window);% longitude edge crossing has to be addressed
    end
    %% seam crossing?
    seam=false;
    if ~full_globe.x
        if (Wlin.west==1  && Wlin.east==size(IN.LON,2)) % ie not full globe but both seam ends are within desired window
            seam=true; % piece crosses long seam
            [OUT,window]=SeamCross(IN,window);
        else % desired piece is within global fields, not need for stitching
            OUT=AllGood(IN,Wlin);
        end
    end
    %% append params
    OUT.window				 =window;
    OUT.params.full_globe =full_globe;
    OUT.params.seam       =seam;
end
function [OUT]=AppendIfFullZonal(IN,window)
    %% append 1/10 of map to include eddies on seam
    % S04_track_eddies is able to avoid counting 1 eddy twice
    ss=window.limits.south;
    nn=window.limits.north;
    %% init
    OUT.grids.LON=IN.LON(ss:nn,:);
    OUT.grids.LAT=IN.LAT(ss:nn,:);
    OUT.grids.SSH=IN.SSH(ss:nn,:);
    %% append
    xadd=round(window.size.X/10);
    OUT.grids.LON=OUT.grids.LON(:,[1:end, 1:xadd]);
    OUT.grids.LAT=OUT.grids.LAT(:,[1:end, 1:xadd]);
    OUT.grids.SSH=OUT.grids.SSH(:,[1:end, 1:xadd]);
end
function [OUT,window]=SeamCross(IN,window)
    Wflag=window.flag;
    Wlin=window.limits;
    %% find new west and east
    easti =find(sum(double(Wflag))==0,1,'first');
    westi =find(sum(double(Wflag))==0,1,'last');
    southi=Wlin.south;
    northi=Wlin.north;
    %% reset east and west
    window.limits.west=westi;
    window.limits.east=easti;
    %% stitch 2 pieces 2g4
    %OUT.window.flag=[Wflag(southi:northi,westi:end) Wflag(southi:northi,1:easti)];
    OUT.grids.LON =[IN.LON(southi:northi,westi:end) IN.LON(southi:northi,1:easti)];
    OUT.grids.LAT =[IN.LAT(southi:northi,westi:end) IN.LAT(southi:northi,1:easti)];
    OUT.grids.SSH =[IN.SSH(southi:northi,westi:end) IN.SSH(southi:northi,1:easti)];
end
function OUT=AllGood(IN,Wlin)
    %% cut piece
    OUT.window.limits=Wlin;
    %OUT.window.flag =Wflag(Wlin.south:Wlin.north,Wlin.west:Wlin.east);
    OUT.grids.LON =IN.LON(Wlin.south:Wlin.north,Wlin.west:Wlin.east);
    OUT.grids.LAT =IN.LAT(Wlin.south:Wlin.north,Wlin.west:Wlin.east);
    OUT.grids.SSH =IN.SSH(Wlin.south:Wlin.north,Wlin.west:Wlin.east);
end
function out=nanLand(in,fac)
    %% nan and SI
    out=in / fac;
    out(out==0)=nan;
end
function [DY,DX]=DYDX(LAT,LON)
    %% grid increment sizes
    DY=deg2rad(abs(diff(double(LAT),1,1)))*earthRadius;
    DX=deg2rad(abs(diff(double(LON),1,2)))*earthRadius.*cosd(LAT(:,1:end-1));
    %% append one line/row to have identical size as other fields
    DY=DY([1:end end],:);
    DX=DX(:,[1:end end]);
    %% correct 360° crossings
    seamcrossflag=DX>100*median(DX(:));
    DX(seamcrossflag)=abs(DX(seamcrossflag) - 2*pi*earthRadius.*cosd(LAT(seamcrossflag)));
end
function file=SampleFile(DD)
    dir_in =DD.path.raw;
    pattern_in=DD.map.pattern.in;
    sample_time=DD.time.from.str;
    file.in=[dir_in.name, strrep(pattern_in, 'yyyymmdd',sample_time)];
    if ~exist(file.in,'file')
        error([file.in,' doesnt exist! choose other start date!'])
    end
end
function [file,exists]=GetCurrentFile(TT,DD) 
    exists.out=false;
   file.in=TT.files;  
   timestr=datestr(TT.daynums,'yyyymmdd');
    %% set up output file
    path=DD.path.cuts.name;
    geo=DD.map.geo;
    file.out=strrep(DD.pattern.in	,'SSSS',sprintf('%04d',geo.south) );
    file.out=strrep(file.out, 'NNNN',sprintf('%04d',geo.north) );
    file.out=strrep(file.out, 'WWWW',sprintf('%04d',geo.west) );
    file.out=strrep(file.out, 'EEEE',sprintf('%04d',geo.east) );
    file.out=[path, strrep(file.out, 'yyyymmdd',timestr)];
    if exist(file.out,'file'), disp([file.out ' exists']); exists.out=true; end
end
function WriteFileOut(file,CUT) %#ok<INUSD>
    save(file,'-struct','CUT')
end

