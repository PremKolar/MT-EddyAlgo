%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 10-Oct-2013 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DD=initialise(toCheck,parentFunc)
	%%
	preInits
	%% get user input
	DD = get_input;
	%% check whether info file exists already
	DDcheck=[DD.path.root, 'DD.mat'];
	
	%% if DD.mat exists, rehash or keep initial mapDims
	
	if exist(DDcheck,'file')
		if DD.switches.rehashMapDims
			DD=catstruct(load(DDcheck),DD);
		else
			DD=catstruct(load(DDcheck),rmfield(DD,'map'));
		end
	end
	%%
	if ~isempty(toCheck)
		DD=ini(DD,toCheck);
	end
	
	%% in case DD was deleted after S00 was executed rehash window info from cut file
	if ~isempty(DD.path.cuts.files)
		[DD.map.window]=GetWin(DD);
	end
	%% load workers
	DD.threads.num=init_threads(DD.threads.num);
	if DD.threads.num>DD.time.span/DD.time.delta_t
		error(toomanythreads,'too many threads for not enough timesteps!!!')
	end
	%% monitoring stuff
	monitorStuff
	%% db stuff
	dbStuff
	%----------------------------------------------------------------------
	%----------------------------------------------------------------------
	function preInits		
		addpath(genpath('./'));  %#ok<*MCAP>
		%         warning on backtrace;
		warning('off','SNCTOOLS:nc_getall:dangerous');
		rehash; clc; close all;
		format shortg;
		dbstop if error;
	end
	%----------------------------------------------------------------------
	function monitorStuff
		DD.monitor.tic=tic;
		if nargin>=2
			DD.monitor.rootFunc=functions(eval(['@' parentFunc]));
		end
		%     dispmem;
	end
	%----------------------------------------------------------------------
	function dbStuff       
		echo off all; diary off;
		if DD.debugmode
			%echo on all
% 			diary on;
			dbstop if error;
		else
			
			for tt=1:DD.threads.num
				commFile=sprintf('./.comm%03d.mat',tt);
				comm=matfile(commFile,'writable',true);
				comm.printstack(1,1)={['thread '  num2str(tt) ]};
			end
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DD=ini(DD,toCheck)
	%% check data for consistency
	DD = check_data(DD,toCheck);
	%% distro thread limits
	maxthreads=feature('numCores');
	threads=min([maxthreads, DD.threads.num]);
	DD.threads.lims = thread_distro(threads,DD.checks.passedTotal);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TT=initChecks(DD,toCheck)
	throwNoDataErr
    %% get filenames
	TT = DD.time;
	TT.existant.filesall=extractfield(DD.path.(toCheck).files,'name');
	%% cat numbers in filenames only for speed
	dateOnly.s=regexp(cat(2,TT.existant.filesall{:}),'\d{8}','match');
	dateOnly.n=cellfun(@(c) datenum(c,'yyyymmdd'),dateOnly.s);
	TT.existant.fcats=cell2mat(cellfun(@(c) ['-' c],dateOnly.s,'uniformoutput',false));
	%% correct start/end date if not exactly on existing time step or not within range
	[~,ii]=min(abs(dateOnly.n-TT.from.num));
	offset=   dateOnly.n(ii) - TT.from.num ;
	TT.from.num=  TT.from.num + offset;
	TT.till.num= TT.till.num + offset; % shift till also to keep const time span
	TT.from.str= datestr(TT.from.num,'yyyymmdd');
	%% correct-shift end of date range also
	[~,ii]=min(abs(dateOnly.n-TT.till.num));
	offset=   dateOnly.n(ii) - TT.till.num ;
	TT.till.num= TT.till.num + offset;
	TT.till.str= datestr(TT.till.num,'yyyymmdd');
	TT.span = TT.till.num - TT.from.num + 1;
	disp(['corrected date range to ' datestr(TT.from.num) ' - ' datestr(TT.till.num)])
	%% init time-steps vector
	TT.timesteps.n = TT.from.num:TT.delta_t:TT.till.num;
	TT.timesteps.s =datestr(TT.timesteps.n,'yyyymmdd');
	%% init new delta t and success vector
	TT.passed = false(numel(TT.timesteps.n),1);
	TT.del_t = nan(size(TT.passed));
	%-----------------------------------------------------------------------
	function throwNoDataErr
		if isempty(DD.path.(toCheck).files)
			error('err:noFiles',['no files found in ' DD.path.(toCheck).name])
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD] = check_data(DD,toCheck)
	%% init
	DD.time=initChecks(DD,toCheck);
	%% check for each needed file
	DD.time.passed=checkForFiles(DD.time);
	checks.del_t =newDt(DD.time);% 'backwards' del_t's
	%% append info
	checks.passedTotal = sum(DD.time.passed);
	checks.passed(checks.passedTotal)=struct;
	temp=num2cell(DD.time.timesteps.n(DD.time.passed))';
	[checks.passed.daynums] =     deal(temp{:});
	%% find corresponding filenames
	checks.passed=getFnames(DD,checks,toCheck);
	%% disp found files
	filedisps(checks);
	%% append
	DD.checks=checks;
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
function del_t=newDt(TT)
	del_t=nan(TT.span,1);
	tempdelt=TT.delta_t;
	for tt = 2:numel(TT.passed);
		if ~TT.passed(tt)
			del_t(tt)=nan;
			tempdelt=tempdelt+TT.delta_t;
		else
			del_t(tt)=tempdelt;
			tempdelt=TT.delta_t;
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
	timestr=cellfun(@(x) datestr(x,'yyyymmdd'),{checks.passed.daynums}','uniformoutput',false); % only passed ones
	for cc=1:numel(timestr)
		ts=timestr{cc};
		if strcmp(toCheck,'raw') % raw filenames relevant
			passed(cc).filenames=[path.name, strrep(DD.map.in.fname, 'yyyymmdd',ts)];
			passed(cc).protofilenames=[];
		else % build new filenames
			geo=DD.map.out;
			file.out=strrep(strrep(pattern, 'yyyymmdd',ts),'CUT',DD.pattern.prefix.(toCheck));
			passed(cc).filenames=[ NSWE2nums(path.name,file.out,geo,ts)  ];
			ii=strfind(passed(cc).filenames,'_');
			passed(cc).protofilenames=passed(cc).filenames(ii:end);
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
