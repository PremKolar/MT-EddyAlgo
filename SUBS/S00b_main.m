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
            disp('exists');continue
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
    %% get data
    for kk={'lat','lon','ssh'}
        keys.(kk{1})=DD.map.in.keys.(kk{1});
    end
    [raw_fields,unreadable]=GetFields(file.in,keys);
    if unreadable.is, CUT.crpt=true; return; end
    %% cut
    [CUT]=cutSlice(raw_fields,DD.map.window);
    %% nan out land and make SI
    CUT.fields.ssh=nanLand(CUT.fields.ssh,DD.parameters.ssh_unitFactor);
    %% get distance fields
    [CUT.fields.dy,CUT.fields.dx]=dydx(CUT.fields.lat,CUT.fields.lon);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=cutSlice(raw,win)
    cutOut=@(raw,idx) reshape(raw(idx),size(idx));
    %% cut piece
    fields=fieldnames(raw);
    for ff=1:numel(fields); field=fields{ff};
        out.fields.(field) = cutOut(raw.(field),win.idx);
    end
    %% append params
    out.window = rmfield(win,{'flag','iy','ix','idx'});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=nanLand(in,fac)
    %% nan and SI
    out=in / fac;
    out(out==0)=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dy,dx]=dydx(lat,lon)
    betweenNodesX = @(lalo) (lalo(:,2:end) + lalo(:,1:end-1))/2;
    betweenNodesY = @(lalo) (lalo(2:end,:) + lalo(1:end-1,:))/2;
    copyBndryX    = @(X) X(:,[1 1:end end]);
    copyBndryY    = @(Y) Y([1 1:end end],:);
    deg2m         = @(degs) deg2km(degs) * 1e3; 
    %% y
    dy=deg2m(abs(diff(lat,1,1)));
    %% x
    dlon=abs(diff(lon,1,2));
    dlon(dlon>180) = abs(dlon(dlon>180) - 360);
    dx=deg2m(dlon) .* cosd(betweenNodesX(lat));
    %% mean back to nodes
    dx=copyBndryX(betweenNodesX(dx));    
    dy=copyBndryY(betweenNodesY(dy));
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
    geo=DD.map.in;
    file.out=NSWE2nums(path,DD.pattern.fname,geo,timestr);
    if exist(file.out,'file'), dispM([file.out ' exists']); exists.out=true; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


