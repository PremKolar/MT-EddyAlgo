%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare ssh data
% reads user input from input_vars.m and map_vars.m
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
    [window]=GetWindow3(file,DD.map.in);
    DD.map.window=window;
    save([DD.path.root 'window.mat'],'window');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD)
    %% distro days to threads
    [II]=SetThreadVar(DD);
    %% loop over files
    S00b_main(DD,II);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function file=SampleFile(DD)
    dir_in =DD.path.raw;
    pattern_in=DD.map.in.fname;
    readable=false;
    sample_time=DD.time.from.str;
    while ~readable
        file=[dir_in.name, strrep(pattern_in, 'yyyymmdd',sample_time)];
        try
            nc_dump(file);
        catch me
            disp(me)
            sample_time=datestr(DD.time.from.num + DD.time.delta_t,'yyyymmdd');
            continue
        end
        readable=true;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



