%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates all contours and saves one file per timestep
function S01a_deleteZonalOverlapInCuts
	%% init
	DD=initialise('cuts');
	%% open pool
	DD.threads.num=init_threads(DD.threads.num);
	%% spmd
	main(DD)
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
function spmd_body(DD)
	%% loop over ssh cuts
	JJ=DD.threads.lims(labindex,1):DD.threads.lims(labindex,2);
	T=disp_progress('init','truncating overlap');
	for jj=JJ
		cutOff(jj,DD);
		T=disp_progress('disp',T,numel(JJ),100);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cutOff(jj,dd)
	%% get data
	[fname,CUT]=getData(jj,dd);
	FN=fieldnames(CUT.grids);
	if ~strcmp(CUT.window.type,'globe'),disp('no need!'); return; end
	%% get actual xsize
	X.right=CUT.window.size.X;
	for fn=FN';fn=fn{1};
		%% cut off
		CUT.grids.(fn)=CUT.grids.(fn)(:,1:X.right);
	end
	CUT.params.truncated=true;
	%% save CUT
	save(fname.full,'-struct','CUT');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fname,CUT]=getData(jj,dd)
	%% load cut
	fname=get_file(jj,dd);
	CUT=load(fname.full);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function file=get_file(jj,dd)
	%% load cut
	file.file=dd.path.cuts.files(jj).name;
	file.base=dd.path.cuts.name;
	file.full =[file.base, file.file ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
