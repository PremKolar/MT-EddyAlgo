%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Sep-2014 04:00:00
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare ssh data
% reads user input from input_vars.m and map_vars.m
function S00A_rho2ssh
     %% init dependencies
    addpath(genpath('./'))
    %% get user input
    DD = initialise('raw',mfilename);      
   
    %%
    DD.path.rawFromRho.name = [DD.path.root 'ssh/'];    
    mkdirp(DD.path.rawFromRho.name);    
    MAP.lat=nc_varget(DD.map.in.LatLonDepthFile,DD.map.in.keys.lat);
    MAP.lon=nc_varget(DD.map.in.LatLonDepthFile,DD.map.in.keys.lat);
    MAP.depth=nc_varget(DD.map.in.LatLonDepthFile,DD.map.in.keys.lat);
    main(DD,MAP)
    %% save info
    conclude(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,MAP)
    if DD.debugmode
        spmd_body(DD,MAP);
    else
        spmd(DD.threads.num)
            spmd_body(DD,MAP);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD,MAP)
    %% distro days to threads
    [II]=SetThreadVar(DD);
    %% loop over files
    S00A_main(DD,II,MAP);
end
function S00A_main(DD,II,MAP)
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
        [SSH]=getSSH(file,DD,MAP);
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
function [ssh]=getSSH(file,DD,MAP)
    %     addpath(genpath('./'));
    %% get data
    for kk={'lat','lon','ssh','rho'}
        keys.(kk{1})=DD.map.in.keys.(kk{1});
    end
    
    ssh.(keys.lat) = MAP.(keys.lat);
    ssh.(keys.lon) = MAP.(keys.lon);
    
    ssh.sourceInfo=nc_info(file.in);
    disp(ssh.sourceInfo);
    rho=nc_varget(file.in,keys.rho, [DD.parameters.zLevel-1 0 0],[1 inf inf]);
   
    
  
end

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
function [file,exists]=GetCurrentFile(TT,DD)
    exists.out=false;
    file.in=TT.files;
    timestr=datestr(TT.daynums,'yyyymmdd');
    %% set up output file
    PATH=DD.path.rawFromRho.name;
    geo=DD.map.out;
    file.out=NSWE2nums(PATH,DD.pattern.fname,geo,timestr);
    if exist(file.out,'file'), dispM([file.out ' exists']); exists.out=true; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


