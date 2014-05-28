function DD=initialise(toCheck)
    addpath(genpath('./'));
    clc;
    close all;
    INPUT = get_input;
    DDcheck=[INPUT.path.root, 'DD.mat'];
    if ~exist('toCheck','var')
        toCheck=false;
    end
    if exist(DDcheck,'file')
        DD=catstruct(load(DDcheck),ini(toCheck,INPUT));
    else
        DD=ini(toCheck,INPUT);
    end
    dbstop if error
    rehash
    format shortg
    
    if ~isfield(DD,'map') && ~isempty(DD.path.cuts.files)
        DD.map=map_vars;
        [DD.map.window]=GetWindow(DD);
    end
end


function INPUT=ini(toCheck,INPUT)
    %% check data for consistency
    if toCheck
        INPUT.checks = check_data(INPUT,toCheck);
        %% distro thread limits
        maxthreads=feature('numCores');
        threads=min([maxthreads, INPUT.threads.num]);
        INPUT.threads.lims = thread_distro(maxthreads,INPUT.checks.passed.total);
    end
end

function checks = check_data(DD,toCheck)
    %% init
    TT = DD.time;
    passed = false(TT.span,1);
    all_time_steps = TT.from.num:TT.delta_t:TT.till.num;
    all_time_steps_str =datestr(all_time_steps,'yyyymmdd');
    all_files=extractfield(DD.path.(toCheck).files,'name');
    %% cat numbers in filenames only for speed
    all_files_cat=cell2mat(regexp(cat(2,all_files{:}),'[0-9]','match'));
    %% init new delta t
    del_t = ones(size(passed))*DD.time.delta_t; del_t(1)=nan;
    %% check for each needed file
    pp = 0;
    T=disp_progress('init','checking data');
    for tt = all_time_steps;
        T=disp_progress('disp',T,numel(all_time_steps),5);
        if (pp>0 && ~passed(pp) && pp<numel(passed))
            del_t(pp+1)=del_t(pp)+ DD.time.delta_t;  % cumsum time steps for missing files
            del_t(pp)=nan;  % nan out del_t for inexistent files
        end
        pp=pp+1;
        current_day = all_time_steps_str(pp,:); 
        % 		current_day = datestr(tt,'yyyymmdd');
        if isempty(findstr(current_day,all_files_cat))           
            continue
        end
        passed(pp)=true;
    end
    
    %% create new del_t time vector in accordance with missing files
    checks.del_t = del_t; % 'backwards' del_t
    
    %% append info
    checks.passed.daynums = all_time_steps(passed)';
    % 	checks.passed.flags = passed;
    checks.passed.total = sum(passed);
    
    % 	%% geo info
    % 	for file=all_files
    % 		w=regexpi(file{1},'.[0-9][0-9][0-9]w')	;
    % 		e=regexpi(file{1},'.[0-9][0-9][0-9]e');
    % 		s=regexpi(file{1},'.[0-9][0-9][0-9]s')	;
    % 		n=regexpi(file{1},'.[0-9][0-9][0-9]n');
    % 		checks.west=str2double(file{1}(w:w+3));
    % 		checks.east=str2double(file{1}(e:e+3));
    % 		checks.south=str2double(file{1}(s:s+3));
    % 		checks.north=str2double(file{1}(n:n+3));
    % 	end
end


function [window]=GetWindow(DD)
    smplFile=[DD.path.cuts.name DD.path.cuts.files(1).name];
    load(smplFile,'window');
    window.flag=[];
end
