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
    [MadFile,CUT]=madsData('../psvar.cdf');
    %% get geo stuff
    [DD,CUT]=geostuff(CUT,DD);
    %% thread distro
    DD.threads.lims=thread_distro(DD.threads.num,numel(CUT.TIME));
    %% create out dir
    [~,~]=mkdir(DD.path.cuts.name);
    %% spmd
    main(DD,MadFile,CUT)
end

function main(DD,MadFile,CUT)
    if DD.debugmode
        spmdBlock(DD,MadFile,CUT)
    else
        spmd(DD.threads.num)
            spmdBlock(DD,MadFile,CUT)
        end
    end
end


function spmdBlock(DD,MadFile,CUT)
    CC=(DD.threads.lims(labindex,1):DD.threads.lims(labindex,2));
    E325=nc_varget(MadFile,'E325');
    %% loop over files
    [T]=disp_progress('init','preparing raw data');
    for cc=CC
        [T]=disp_progress('calc',T,numel(CC),4242);
        operateDay(squeeze(E325(cc,:,:)),CUT,DD,cc);
    end
end
function [MadFile,CUT]=madsData(MadFile)
    try
        nc_info(MadFile);
    catch
        error(['expecting  ' MadFile])
    end
    CUT.TIME=nc_varget(MadFile,'TIME')+datenum('1900','yyyy');
    
    CUT.XT_bnds=nc_varget(MadFile,'XT_bnds');
    CUT.YT_bnds=nc_varget(MadFile,'YT_bnds');
    CUT.XT=nc_varget(MadFile,'XT');
    CUT.YT=nc_varget(MadFile,'YT');
end
function [DD,CUT]=geostuff(CUT,DD)
    [CUT.grids.XX,CUT.grids.YY]=meshgrid(CUT.XT,CUT.YT);
    CUT.grids.LAT=rad2deg(CUT.grids.YY./earthRadius);
    CUT.grids.LON=rad2deg(CUT.grids.XX./(cosd(CUT.grids.LAT)*earthRadius));
    [CUT.grids.DY,CUT.grids.DX]=DYDX(CUT.grids.LAT,CUT.grids.LON);
    DD.map.west=min(CUT.grids.LON(:));
    DD.map.east=max(CUT.grids.LON(:));
    DD.map.south=min(CUT.grids.LAT(:));
    DD.map.north=max(CUT.grids.LAT(:));
    [Y,X]=size(CUT.grids.LON);
    DD.map.window.size.X=X;
    DD.map.window.size.Y=Y;
    DD.map.window.limits.west=1;
    DD.map.window.limits.east=X;
    DD.map.window.limits.south=1;
    DD.map.window.limits.north=Y;
end
function operateDay(SSH,CUT,DD,cc)
    %% set up output file
    tt=CUT.TIME(cc);
    timestr=datestr(tt,'yyyymmdd');
    path=DD.path.cuts.name;
    geo=DD.map;
    file.out=strrep(DD.pattern.fname	,'SSSS',sprintf('%04d',round(geo.south)) );
    file.out=strrep(file.out, 'NNNN',sprintf('%04d',round(geo.north) ));
    file.out=strrep(file.out, 'WWWW',sprintf('%04d',round(geo.west) ));
    file.out=strrep(file.out, 'EEEE',sprintf('%04d',round(geo.east)) );
    file.out=[path, strrep(file.out, 'yyyymmdd',timestr)];
    if exist(file.out,'file'), return; end
    %% weird values on borders..
    [Y,X]=size(CUT.grids.LON);
    SSH(1:end,1)=SSH(1:end,2);
    SSH(1:end,X)=SSH(1:end,X-1);
    SSH(1,1:end)=SSH(2,1:end);
    SSH(Y,1:end)=SSH(Y-1,1:end);
    SSH(SSH>10000)=nan;
    SSH(SSH<-10000)=nan;
    CUT.grids.SSH=double(SSH/DD.map.in.SSH_unitFactor);
    %%
    %% append 1/4th zonally to track eddies crossing open bndry
    CUT.params.full_globe.x=true;
    for FN=fieldnames(CUT.grids)';fn=FN{1};
        CUT.grids.(fn)= CUT.grids.(fn)(:,[1:end 1:1/4*X]);
    end
    %% save actual size
    CUT.window.size.X=X;
    CUT.window.size.Y=Y;
    %%
    save(file.out,'-struct','CUT')
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
