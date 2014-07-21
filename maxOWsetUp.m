%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD]=maxOWsetUp
    %% init
    DD=initialise([],mfilename);  
    %% threads
    DD.threads.num=init_threads(DD.threads.num);
    %% find temp and salt files
    [DD.path.TSow]=DataInit(DD);
    %% get window according to user input
    [DD.TSow.window,~]=GetWindow(DD.path.TSow.files(1).salt,DD.map.in,DD.TS.keys);
    %% get z info
    DD.TSow.window=mergeStruct2(DD.TSow.window, GetFields(DD.path.TSow.files(1).salt, cell2struct({DD.TS.keys.depth},'depth')));
    DD.TSow.window.size.Z=numel(DD.TSow.window.depth);
    %% distro time steps to threads
    DD.TSow.lims=thread_distro(DD.threads.num,numel(DD.path.TSow.files));
    
    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out]=DataInit(DD)
    %% find the temp and salt files
DirIn=dir([DD.path.full3d.name '*.nc'])   ;
tt=0;ss=0;
    for kk=1:numel(DirIn);
        if ~isempty(strfind(upper(DirIn(kk).name),DD.TS.keys.salt))
            ss=ss+1;
            out.files(ss).salt=[DD.path.full3d.name DirIn(kk).name];
        end
        if ~isempty(strfind(upper(DirIn(kk).name),DD.TS.keys.temp))
            tt=tt+1;
            out.files(tt).temp=[DD.path.full3d.name DirIn(kk).name];
        end
    end
    
     out.fnum=numel(out.files);
    out.dir   = [DD.path.Rossby.name];  
out.dailyBaseName   = [ out.dir 'OW_']; 
    
end
