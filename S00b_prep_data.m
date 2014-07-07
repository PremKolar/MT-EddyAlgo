%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare ssh data
% reads user input from input_vars.m and map_vars.m
% This is the only file that would have to adapted to the strcture of the
% input ssh data. This step is not officially part of the program. Use a
% copy of this and adapt to your data, so that S01 gets the required input
% structure.
function S00b_prep_data
    %% set up
    [DD]=set_up;
    %% spmd
    main(DD)
    %% save info
    conclude(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD)
    if DD.debugmode
        spmd_body(DD);
    else
        spmd(DD.threads.num)
            spmd_body(DD);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD]=set_up
    %% init dependencies
    addpath(genpath('./'))
    %% get user input
    DD = initialise('raw',mfilename);
    %% get sample window
    file=SampleFile(DD);
    [DD.map.window]=GetWindow2(file,DD.map.in);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD)
    %% distro days to threads
    [II]=SetThreadVar(DD);
    %% loop over files
    [T]=disp_progress('init','preparing raw data');
    for cc=1:numel(II);
        %%
        [T]=disp_progress('calc',T,numel(II),100);
        %% get data
        file=GetCurrentFile(II(cc),DD)  ;
        %% cut data
        [CUT]=CutMap(file,DD);  % if isempty(CUT); return; end
        %% filter
        CUT.grids.sshS=lowPassSSH(CUT.grids,DD.parameters.Gausswidth);
        %% write data
        WriteFileOut(file.out,CUT);
    end
    disp_progress('conclude');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CUT]=CutMap(file,DD)
    addpath(genpath('./'));
    %% get data
    for kk={'lat','lon','ssh'}
        keys.(kk{1})=DD.map.in.keys.(kk{1});
    end
    [raw_fields,unreadable]=GetFields(file.in,keys);
    if unreadable.is, CUT=[]; return; end
    %% cut
    [CUT]=ZonalAppend(raw_fields,DD.map.window);
    %% nan out land and make SI
    CUT.grids.ssh=nanLand(CUT.grids.ssh,DD.parameters.ssh_unitFactor);
    %% get distance fields
    [CUT.grids.DY,CUT.grids.DX]=DYDX(CUT.grids.lat,CUT.grids.lon);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OUT=ZonalAppend(raw,window)    
    OUT=buildGrids;
    %----------------------------------------------------------------------
    function out=buildGrids
        xlin=drop_2d_to_1d(window.iy,window.ix,window.fullsize(1));
        %% cut piece
        fields=fieldnames(raw);
        for ff=1:numel(fields); field=fields{ff};
            out.grids.(field) = raw.(field)(xlin);
        end
        %% append params
        out.window = rmfield(window,{'flag','iy','ix'});       
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=lowPassSSH(grids,gaussWidth)

gw.x=gaussWidth./grids.DX;
    gw.y=gaussWidth./grids.DY;
    
    % FOR NOW    
    gwx=ceil(mean(gw.x(:)));
    gwy=ceil(mean(gw.y(:)));
    
%     a=smooth2gauss(grids.ssh,gwx,gwy)
    
    
   out = smooth2gauss(grids.ssh,gwx,gwy);
%     m = m - nanmin(m(:));
%     m = m/nanmax(m(:));
%     out = m*(nanmax(grids.ssh(:)) - nanmin(grids.ssh(:))) + nanmin(grids.ssh(:));
%     mm=smooth2a(grids.ssh,3,3)
%     
%      figure(5);
%      ppcolor(mm)
% %      caxis([-.3 .3])
%      figure(1);
%      ppcolor(grids.ssh)
% %      caxis([-.3 .3])
% figure(2);
%     ppcolor(m)
% % caxis([-.3 .3])
% 
% figure(3);
% contour(m,[0 0])
% hold on
% contour(grids.ssh,[0 0],'color',[1 0 0])
% 
% figure(4);
% subplot(1,3,1)
% contourf(m,120)
% subplot(1,3,2)
% contourf(grids.ssh,120)
% subplot(1,3,3)
% contourf(mm,120)
% 

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
function file=SampleFile(DD)
    dir_in =DD.path.raw;
    pattern_in=DD.map.in.fname;
    sample_time=DD.time.from.str;
    file=[dir_in.name, strrep(pattern_in, 'yyyymmdd',sample_time)];
    if ~exist(file,'file')
        error([file,' doesnt exist! choose other start date!'])
    end
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


