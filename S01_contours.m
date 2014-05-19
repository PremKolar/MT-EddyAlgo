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
%     spmd(DD.threads.num)
        spmd_body(DD);
%     end
    %% save info
    save_info(DD)
end
function spmd_body(DD)
    DD.id=labindex;
    %% loop over ssh cuts
    JJ=DD.threads.lims(DD.id,1):DD.threads.lims(DD.id,2);
    II=struct;
    for jj=JJ
        %% contours
        II=get_contours(jj,DD,JJ,II);
    end
end
function II=get_contours(jj,dd,JJ,II)
    %% check
    outFile=[dd.path.conts.name regexprep(dd.path.cuts.files(jj).name, 'CUT', 'CONT')];
    if exist(outFile,'file')
        return
    end
    %%
    [II,CONT]=init_get_contours(jj,dd,JJ,outFile);
    %% loop over levels
    for level=II.levels
        II.T=disp_progress('disp',II.T,numel(II.levels),3);
        CONT.all=[CONT.all; contourc(II.grids.SSH,[level level])'];
    end
    %% save data
    save(CONT.filename,'-struct','CONT');
end
function [II,CONT]=init_get_contours(jj,dd,JJ,filename)
    %% load cut
    II.file=get_file(jj,dd);
    II.grids=getfield(load(II.file.full),'grids');
    %% calc contours
    disp('calculating contours... takes long time!')
    CONT.all=[]; % init
    %% create level vector at chosen interval
    II.levels=nanmin(II.grids.SSH(:))-0.1:dd.contour.step:nanmax(II.grids.SSH(:))+0.1;
    II.T=disp_progress('init',['contours of day: ',[sprintf('%03i',jj+1-dd.threads.lims(dd.id,1)),'/',sprintf('%03i',numel(JJ))]]);
    II.T.uplevel.togo=numel(JJ)-(jj-JJ(1));
    II.T.uplevel.perDone=[ num2str(round((jj+1-JJ(1))/numel(JJ)*100)),' %'];
    II.T.uplevel.l=(jj+1-(JJ(1)));
    II.T.uplevel.L=numel(JJ);
    %% add info
    CONT.filename=filename;
end
function file=get_file(jj,dd)
    %% load cut
    file.file=dd.path.cuts.files(jj).name;
    file.base=dd.path.cuts.name;
    file.full =[file.base, file.file ];
end

