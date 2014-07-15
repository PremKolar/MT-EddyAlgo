%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD]=maxOWsetUp
	addpath(genpath('../'))
	addpath(genpath('../SUBS/'))	
	%% init
	DD=initialise([],mfilename);
	%% check if exists already
	[DD.path.Rossby.NCfile] = initNC(DD);
	%% threads
	DD.threads.num=init_threads(DD.threads.num);
	%% find temp and salt files
	[DD.path.TempSalt]=tempsalt(DD);
	%% get window according to user input
	[DD.TS.window,~]=GetWindow(DD.path.TempSalt.salt(1),DD.map.in,DD.TS.keys);
	%% distro X lims to chunks
	DD.RossbyStuff.lims.data=limsdata(DD.parameters.RossbySplits,DD.TS.window);
	%% distro chunks to threads
	DD.RossbyStuff.lims.loop=thread_distro(DD.threads.num,DD.parameters.RossbySplits);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lims=limsdata(splits,window)
	%% set dimension for splitting (files dont fit in memory)
	X=window.size.X;
	%% distro X lims to chunks
	lims=thread_distro(splits,X) + window.limits.west-1;
	%% in case window crosses zonal bndry
	lims(lims>window.fullsize(2)) = lims(lims>window.fullsize(2)) - window.fullsize(2);
	%% in case one chunk crosses zonal bndry
	td=lims(:,2)-lims(:,1) < 0; % find chunk
	lims(td,1)=0; % let it start at 0
	lims(find(td)-1,2)=window.fullsize(2)-1; % let the one before finish at end(X)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outfilename] = initNC(DD)
	outfilename=[DD.path.Rossby.name, 'OW.nc'];
	NCoverwriteornot(outfilename);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [file]=tempsalt(DD)
	%% find the temp and salt files
	tt=0;ss=0;
	for kk=1:numel(DD.path.TempSalt.files);
		if ~isempty(strfind(upper(DD.path.TempSalt.files(kk).name),'SALT'))
			ss=ss+1;
			file(ss).salt=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name];
		end
		if ~isempty(strfind(upper(DD.path.TempSalt.files(kk).name),'TEMP'))
			tt=tt+1;
			file(tt).temp=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name];
		end
	end
end
