function DD=initialise(toCheck)
	addpath(genpath('./'));
	clc;
	close all;
	temp = get_input;
	DDcheck=[temp.path.root, 'DD.mat'];
	if ~exist('toCheck','var')
		toCheck=false;
	end
	
	if exist(DDcheck,'file')
		DD=catstruct(load(DDcheck),ini(toCheck));
	else
		DD=ini(toCheck);
	end	
	dbstop if error
	rehash	
end

function out=ini(toCheck)
	%% read input file
	out = get_input;
	%% check data for consistency
	if toCheck
		out.checks = check_data(out,toCheck);	
	%% distro thread limits
	out.threads.lims = thread_distro(out.threads.num,out.checks.passed.total);
	end
end

function checks = check_data(DD,toCheck)
	%% init
	TT = DD.time;
	passed = false(TT.span,1);
	all_time_steps = TT.from.num:TT.delta_t:TT.till.num;
	all_files=struct2cell(DD.path.(toCheck).files);
	all_files=all_files(1,:); % cell with filenames
	all_passed = false(numel(all_files),1);
	%% check for each needed file
	pp = 0;
	%% init new delta t
	del_t = ones(size(passed))*DD.time.delta_t; del_t(1)=nan;
	T=disp_progress('init','checking data');
	for tt = all_time_steps;
	T=disp_progress('disp',T,numel(all_time_steps),10);		
	if (pp>0 && ~passed(pp) && pp<numel(passed))
			del_t(pp+1)=del_t(pp)+ DD.time.delta_t;  % cumsum time steps for missing files
			del_t(pp)=nan; % nan out del_t for inexistent files			
	end	
		pp=pp+1;
		current_day = datestr(tt,'yyyymmdd');
		for aa=1:numel(all_passed)
			if isempty(findstr(current_day,all_files{aa}))
				continue
			end
			all_passed(aa)=true;
			passed(pp)=true;
		end
	end
	%% create new del_t time vector in accordance with missing files
	checks.del_t = del_t; % 'backwards' del_t
	
	%% append info
	checks.passed.daynums = all_time_steps(passed)';
	checks.passed.flags = passed;
	checks.passed.total = sum(passed);
	checks.files=all_files(all_passed);
	
	%% geo info
	for file=all_files
		w=regexpi(file{1},'.[0-9][0-9][0-9]w')	;
		e=regexpi(file{1},'.[0-9][0-9][0-9]e');
		s=regexpi(file{1},'.[0-9][0-9][0-9]s')	;
		n=regexpi(file{1},'.[0-9][0-9][0-9]n');
		checks.west=str2double(file{1}(w:w+3));
		checks.east=str2double(file{1}(e:e+3));
		checks.south=str2double(file{1}(s:s+3));
		checks.north=str2double(file{1}(n:n+3));
	end	
end
