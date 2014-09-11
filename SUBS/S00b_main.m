%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Sep-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub to ../S00b_prep..
function S00b_main(DD,II)
    [T]=disp_progress('init','preparing raw data');
    for cc=1:numel(II);
        [T]=disp_progress('calc',T,numel(II),100);
        %% get data
        [file,exists]=GetCurrentFile(II(cc),DD)  ;
        %% skip if exists and ~overwrite switch
        if exists.out && ~DD.overwrite;
            disp('exists');return
        end
        %% cut data
        [CUT]=CutMap(file,DD);
        %% save empty corrupt files too
        if isfield(CUT,'crpt');
            [d,f,x] = fileparts(file.out ) ;
            file.out = fullfile(d,['CORRUPT-' f x]);
        end
        %% write data
        WriteFileOut(file.out,CUT);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CUT]=CutMap(file,DD)
    %     addpath(genpath('./'));
    %% get data
    for kk={'lat','lon','ssh'}
        keys.(kk{1})=DD.map.in.keys.(kk{1});
    end
    [raw_fields,unreadable]=GetFields(file.in,keys);
    if unreadable.is, CUT.crpt=true; return; end
    %% cut
    [CUT]=ZonalAppend(raw_fields,DD.map.window);
    %% nan out land and make SI
    CUT.grids.ssh=nanLand(CUT.grids.ssh,DD.parameters.ssh_unitFactor);
    %% get distance fields
    [CUT.grids.DY,CUT.grids.DX]=DYDX(CUT.grids.lat,CUT.grids.lon);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=ZonalAppend(raw,window)
    xlin=drop_2d_to_1d(window.iy,window.ix,window.fullsize(1));
    %% cut piece
    fields=fieldnames(raw);
    for ff=1:numel(fields); field=fields{ff};
        out.grids.(field) = raw.(field)(xlin);
    end
    %% append params
    out.window = rmfield(window,{'flag','iy','ix'});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=nanLand(in,fac)
    %% nan and SI
    out=in / fac;
    out(out==0)=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DY,DX]=DYDX(lat,lon)
    %% grid increment sizes
    DY=deg2rad(abs(diff(double(lat),1,1)))*earthRadius;
    DX=deg2rad(abs(diff(double(lon),1,2)))*earthRadius.*cosd(lat(:,1:end-1));
    %% append one line/row to have identical size as other fields
    DY=DY([1:end end],:);
    DX=DX(:,[1:end end]);
    %% correct 360° crossings
    seamcrossflag=DX>100*median(DX(:));
    DX(seamcrossflag)=abs(DX(seamcrossflag) - 2*pi*earthRadius.*cosd(lat(seamcrossflag)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========================================================================%
function WriteFileOut(file,CUT) %#ok<INUSD>
    save(file,'-struct','CUT')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [file,exists]=GetCurrentFile(TT,DD)
    exists.out=false;
    file.in=TT.files;
    timestr=datestr(TT.daynums,'yyyymmdd');
    %% set up output file
    path=DD.path.cuts.name;
    geo=DD.map.out;
    file.out=NSWE2nums(path,DD.pattern.fname,geo,timestr);
    if exist(file.out,'file'), dispM([file.out ' exists']); exists.out=true; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


