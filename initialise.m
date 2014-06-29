%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 10-Oct-2013 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DD=initialise(toCheck)
	%% very first settings
	addpath(genpath('./'));  %#ok<*MCAP>
	warning on backtrace;
	dbstop if error;
	rehash; clc; close all;
	format shortg;
	%% get user input
	DD = get_input;
	%% check whether info file exists already
	DDcheck=[DD.path.root, 'DD.mat'];
	if ~exist('toCheck','var')
		toCheck=false;
	end
	if toCheck
		DD=ini(DD,toCheck);	
	end
	%% if exist append new info from ini() but keep info not overwritten by ini()
	if exist(DDcheck,'file')
		DD=catstruct(load(DDcheck),DD);
	end
	%% in case DD was deleted after S00 was executed rehash window info from cut file
	if   ~isempty(DD.path.cuts.files)
		[DD.map.window]=GetWin(DD);
	end
	%% load workers
	DD.threads.num=init_threads(DD.threads.num);
	if DD.threads.num>DD.time.span/DD.time.delta_t
		error(toomanythreads,'too many threads for not enough timesteps!!!')
	end
	%% performance stuff
	DD.tic=tic;
	%     dispmem;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DD=ini(DD,toCheck)
	%% check data for consistency
	[DD.checks,abort] = check_data(DD,toCheck);
	if abort,return;end
	%% distro thread limits
	maxthreads=feature('numCores');
	threads=min([maxthreads, DD.threads.num]);
	DD.threads.lims = thread_distro(threads,DD.checks.passedTotal);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TT=initChecks(DD,toCheck)
	TT = DD.time;
	TT.timesteps.n = TT.from.num:TT.delta_t:TT.till.num;
	TT.passed = false(numel(TT.timesteps.n),1);
	TT.timesteps.s =datestr(TT.timesteps.n,'yyyymmdd');
	TT.existant.filesall=extractfield(DD.path.(toCheck).files,'name');
	%% cat numbers in filenames only for speed
	TT.existant.fcats=cell2mat(regexp(cat(2,TT.existant.filesall{:}),'[0-9]','match'));
	%% init new delta t
	TT.del_t = nan(size(TT.passed));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [checks,abort] = check_data(DD,toCheck)
	%% init
	TT=initChecks(DD,toCheck);
	%% check for each needed file
	TT.passed=checkForFiles(TT);
	if ~any(TT.passed)
		checks=nan; abort=true;
		ls(DD.path.(toCheck).name)
		warning(['found no ' toCheck ' files'])
		return
	end
	abort=false;
	%% calc dt's
	TT.del_t=newDt(TT,DD);
	%% append info
	checks.del_t = TT.del_t; % 'backwards' del_t
	checks.passedTotal = sum(TT.passed);
	checks.passed(checks.passedTotal)=struct;
	temp=num2cell(TT.timesteps.n(TT.passed))';
	[checks.passed.daynums] =     deal(temp{:});
	%% find corresponding filenames
	checks.passed=getFnames(DD,checks,toCheck);
	%% disp found files
	filedisps(checks)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filedisps(checks)
	sleep(1)
	disp(['found '])
	for ff=1:numel(checks.passed)
		disp([checks.passed(ff).filenames])
	end
	disp(['total of ' num2str(numel(checks.passed))])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pass=checkForFiles(TT)
	pass=TT.passed;
	for tt = 1:numel(TT.passed);
		if ~isempty(strfind(TT.existant.fcats, TT.timesteps.s(tt,:)))
			pass(tt)=true;
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function del_t=newDt(TT,DD)
	del_t=nan(TT.span,1);
	tempdelt=DD.time.delta_t;
	for tt = 2:numel(TT.passed);
		if ~TT.passed(tt)
			del_t(tt)=nan;
			tempdelt=tempdelt+DD.time.delta_t;
		else
			del_t(tt)=tempdelt;
			tempdelt=DD.time.delta_t;
		end
	end
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window]=GetWin(DD)
	smplFile=[DD.path.cuts.name DD.path.cuts.files(1).name];
	load(smplFile,'window');
	window.flag=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function passed=getFnames(DD,checks,toCheck)
	passed=checks.passed;
	path=DD.path.(toCheck);
	pattern=DD.pattern.fname;
	timestr=cellfun(@(x) datestr(x,'yyyymmdd'),{checks.passed.daynums}','uniformoutput',false);
	cc=0;
	for ts=timestr';cc=cc+1;
		if strcmp(toCheck,'raw')
			passed(cc).filenames=[path.name, strrep(DD.map.in.fname, 'yyyymmdd',ts{1})];
			passed(cc).protofilenames=[];
		else
			temp=[path.name, strrep(pattern, 'yyyymmdd',ts{1})];
			temp=strrep(temp	,'SSSS',sprintf('%04d',DD.map.in.south) );
			temp=strrep(temp	,'NNNN',sprintf('%04d',DD.map.in.north) );
			temp=strrep(temp	,'WWWW',sprintf('%04d',DD.map.in.west) );
			temp=strrep(temp	,'EEEE',sprintf('%04d',DD.map.in.east) );
			passed(cc).filenames=strrep(temp	,'CUT',DD.pattern.prefix.(toCheck));
			ii=strfind(temp,'_');
			passed(cc).protofilenames=temp(ii:end);
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



