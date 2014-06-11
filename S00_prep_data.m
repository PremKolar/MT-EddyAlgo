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
    %% get sample window
     file=SampleFile(DD);
    [DD.map.window]=GetWindow(file,DD.map.in,DD.map.in.keys);
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
        [CUT]=CutMap(file,DD);
        %% write data
        WriteFileOut(file.out,CUT);
    end
end
function [CUT]=CutMap(file,DD)
    addpath(genpath('./'));
    %% get data
    [raw_fields]=GetFields(file.in,DD.map.in.keys);
    %% cut
    [CUT]=ZonalProblem(raw_fields,DD.map.window);
    %% nan out land and make SI
    CUT.grids.ssh=nanLand(CUT.grids.ssh,DD.map.in.ssh_unitFactor);
    %% get distance fields
    [CUT.grids.DY,CUT.grids.DX]=DYDX(CUT.grids.lat,CUT.grids.lon);
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
%=========================================================================%
function WriteFileOut(file,CUT) %#ok<INUSD>
    save(file,'-struct','CUT')
end
function file=SampleFile(DD)
    dir_in =DD.path.raw;
    pattern_in=DD.map.in.fname;
    sample_time=DD.time.from.str;
    file=[dir_in.name, strrep(pattern_in, 'yyyymmdd',sample_time)];
    if ~exist(file,'file')
        error([file,' doesnt exist! choose other start date!'])
    end
end
function [file,exists]=GetCurrentFile(TT,DD) 
    exists.out=false;
   file.in=TT.files;  
   timestr=datestr(TT.daynums,'yyyymmdd');
    %% set up output file
    path=DD.path.cuts.name;
    geo=DD.map.in;
    file.out=strrep(DD.pattern.fname	,'SSSS',sprintf('%04d',geo.south) );
    file.out=strrep(file.out, 'NNNN',sprintf('%04d',geo.north) );
    file.out=strrep(file.out, 'WWWW',sprintf('%04d',geo.west) );
    file.out=strrep(file.out, 'EEEE',sprintf('%04d',geo.east) );
    file.out=[path, strrep(file.out, 'yyyymmdd',timestr)];
    if exist(file.out,'file'), disp([file.out ' exists']); exists.out=true; end
end
