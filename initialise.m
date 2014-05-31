function DD=initialise(toCheck)
    addpath(genpath('./'));
    clc;
    close all;
    INPUT = get_input;
       INPUT.map=map_vars;
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
    
     
    if   ~isempty(DD.path.cuts.files)
        [DD.map.window]=GetWindow(DD);
    end
    
   
    DD.threads.num=init_threads(DD.threads.num   );
end


function INPUT=ini(toCheck,INPUT)
    %% check data for consistency
    if toCheck
        INPUT.checks = check_data(INPUT,toCheck);
        %% distro thread limits
        maxthreads=feature('numCores');
        threads=min([maxthreads, INPUT.threads.num]);
        INPUT.threads.lims = thread_distro(threads,INPUT.checks.passed.total);
    end
end

function checks = check_data(DD,toCheck)
    %% init
    TT = DD.time;
    all_time_steps = TT.from.num:TT.delta_t:TT.till.num;
    passed = false(numel(all_time_steps),1);
    all_time_steps_str =datestr(all_time_steps,'yyyymmdd');
    all_files=extractfield(DD.path.(toCheck).files,'name');
    %% cat numbers in filenames only for speed
    all_files_cat=cell2mat(regexp(cat(2,all_files{:}),'[0-9]','match'));
    %% init new delta t
    del_t = nan(size(passed));
    %% check for each needed file
   for tt = 1:numel(passed);     
         if ~isempty(strfind(all_files_cat, all_time_steps_str(tt,:)))           
             passed(tt)=true;             
        end
   end
   %% calc td's
   tempdelt=DD.time.delta_t;
    for tt = 2:numel(passed);     
        if ~passed(tt)
            del_t(tt)=nan;
            tempdelt=tempdelt+DD.time.delta_t;
        else
            del_t(tt)=tempdelt;
            tempdelt=DD.time.delta_t;
        end
    end     
    %% create new del_t time vector in accordance with missing files
    checks.del_t = del_t; % 'backwards' del_t    
    %% append info
    checks.passed.daynums = all_time_steps(passed)';
    checks.passed.total = sum(passed);   
end

function [window]=GetWindow(DD)
    smplFile=[DD.path.cuts.name DD.path.cuts.files(1).name];
    load(smplFile,'window');
    window.flag=[];    
end
