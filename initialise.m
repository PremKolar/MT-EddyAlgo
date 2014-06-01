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
        INPUT.threads.lims = thread_distro(threads,INPUT.checks.passedTotal);
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
     checks.passedTotal = sum(passed);
    checks.passed(checks.passedTotal)=struct;
    temp=num2cell(all_time_steps(passed))';
    [checks.passed.daynums] =     deal(temp{:});   
    
      %% find corresponding file   
    path=DD.path.(toCheck);
    pattern=DD.pattern.in; 
    timestr=cellfun(@(x) datestr(x,'yyyymmdd'),{checks.passed.daynums}','uniformoutput',false);
    cc=0;
    for ts=timestr';cc=cc+1;
       if strcmp(toCheck,'raw')
         checks.passed(cc).filenames=[path.name, strrep(DD.map.pattern.in, 'yyyymmdd',ts{1})];
       else
        temp=[path.name, strrep(pattern, 'yyyymmdd',ts{1})];
        temp=strrep(temp	,'SSSS',sprintf('%04d',DD.map.geo.south) );
        temp=strrep(temp	,'NNNN',sprintf('%04d',DD.map.geo.north) );
        temp=strrep(temp	,'WWWW',sprintf('%04d',DD.map.geo.west) );
        temp=strrep(temp	,'EEEE',sprintf('%04d',DD.map.geo.east) );
        checks.passed(cc).filenames=strrep(temp	,'CUT',DD.pattern.prefix.(toCheck));  
        ii=strfind(temp,'_');
        checks.passed(cc).protofilenames=temp(ii:end);
       end
    end    
end

function [window]=GetWindow(DD)
    smplFile=[DD.path.cuts.name DD.path.cuts.files(1).name];
    load(smplFile,'window');
    window.flag=[];    
end
