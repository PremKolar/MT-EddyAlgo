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
    [MadFile,MF]=madsData('../psvar.cdf');
    %% get geo stuff
    [DD,MF]=geostuff(MF,DD);
    %% thread distro
    DD.threads.lims=thread_distro(DD.threads.num,numel(MF.TIME));
    %% create out dir
    [~,~]=mkdir(DD.path.cuts.name);
    %% spmd
    main(DD,MadFile,MF)
end

function main(DD,MadFile,MF)
    if DD.debugmode
        spmdBlock(DD,MadFile,MF)
    else
        spmd(DD.threads.num)
            spmdBlock(DD,MadFile,MF)
        end
    end
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
function [MadFile,MF]=madsData(MadFile)
    try
        nc_info(MadFile);
    catch
        error(['expecting  ' MadFile])
    end
    MF.TIME=nc_varget(MadFile,'TIME')+datenum('1900','yyyy');
    MF.XT_bnds=nc_varget(MadFile,'XT_bnds');
    MF.YT_bnds=nc_varget(MadFile,'YT_bnds');
    MF.XT=nc_varget(MadFile,'XT');
    MF.YT=nc_varget(MadFile,'YT');
end
function [DD,MF]=geostuff(MF,DD)
    [MF.grids.XX,MF.grids.YY]=meshgrid(MF.XT,MF.YT);
    MF.grids.LAT=rad2deg(MF.grids.YY./earthRadius);
    MF.grids.LON=rad2deg(MF.grids.XX./(cosd(MF.grids.LAT)*earthRadius));
    [MF.grids.DY,MF.grids.DX]=DYDX(MF.grids.LAT,MF.grids.LON);
    DD.map.west=min(MF.grids.LON(:));
    DD.map.east=max(MF.grids.LON(:));
    DD.map.south=min(MF.grids.LAT(:));
    DD.map.north=max(MF.grids.LAT(:));
    [Y,X]=size(MF.grids.LON);
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
    file.out=strrep(DD.pattern.fname	,'SSSS',sprintf('%04d',round(geo.south)) );
    file.out=strrep(file.out, 'NNNN',sprintf('%04d',round(geo.north) ));
    file.out=strrep(file.out, 'WWWW',sprintf('%04d',round(geo.west) ));
    file.out=strrep(file.out, 'EEEE',sprintf('%04d',round(geo.east)) );
    file.out=[path, strrep(file.out, 'yyyymmdd',timestr)];
    if exist(file.out,'file'), return; end
    %% weird values on borders..
    [Y,X]=size(MF.grids.LON);
    SSH(1:end,1)=SSH(1:end,2);
    SSH(1:end,X)=SSH(1:end,X-1);
    SSH(1,1:end)=SSH(2,1:end);
    SSH(Y,1:end)=SSH(Y-1,1:end);
    SSH(SSH>10000)=nan;
    SSH(SSH<-10000)=nan;
    MF.grids.SSH=double(SSH/DD.map.in.SSH_unitFactor);
    %%
    %% append 1/4th zonally to track eddies crossing open bndry
    MF.params.full_globe.x=true;
    for FN=fieldnames(MF.grids)';fn=FN{1};
        MF.grids.(fn)= MF.grids.(fn)(:,[1:end 1:1/4*X]);
    end
    
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
