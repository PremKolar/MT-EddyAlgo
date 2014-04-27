%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates all contours and saves one file per timestep
function S01_contours
	%% init
	DD=initialise('cuts');
	%% open pool
	init_threads(DD.threads.num)
	%% spmd
	spmd
		spmd_body(DD);
	end
	%% save info
	save_info(DD)
	
end
function spmd_body(DD)
	
	DD.id=labindex;
	%% loop over ssh cuts
	JJ=DD.threads.lims(DD.id,1):DD.threads.lims(DD.id,2);
	for jj=JJ
		%% contours
		get_contours(jj,DD,JJ);
	end
end
function get_contours(jj,dd,JJ)
	%% check
	outFile=[dd.path.conts.name regexprep(dd.path.cuts.files(jj).name, 'CUT', 'CONT')];		
	if exist(outFile,'file')
		return
	end	
	%%
	[II,CONT]=init_get_contours(jj,dd,JJ,outFile);
	%% loop over levels
	for level=II.levels
	dfb
		II.T=disp_progress('disp',II.T,numel(II.levels),5,II.days_prog);
		CONT.all=[CONT.all; contourc(II.grids.SSH,[level level])'];
	end
	%% save data	
	save(CONT.filename,'-struct','CONT');
end
function [II,CONT]=init_get_contours(jj,dd,JJ,filename)
	%% load cut
	II.file=get_file(jj,dd);
	II.grids=getfield(load(II.file.full),'grids');	 %#ok<GFLD>
	%% calc contours
	disp('calculating contours... takes long time!')
	II.days_prog=(jj-JJ(1)+1)/numel(JJ); % %/100 done of days
	CONT.all=[]; % init
	%% create level vector at chosen interval
	II.levels=nanmin(II.grids.SSH(:))-0.1:dd.contour.step:nanmax(II.grids.SSH(:))+0.1;
	II.T=disp_progress('init',['contours of day: ',[sprintf('%03i',jj+1-dd.threads.lims(dd.id,1)),'/',sprintf('%03i',numel(JJ))]]);
	%% add info
	CONT.filename=filename;
	CONT.input=dd;  % add input info
end

function file=get_file(jj,dd)
	%% load cut
	file.file=dd.path.cuts.files(jj).name;
	file.base=dd.path.cuts.name;
	file.full =[file.base, file.file ];
end

