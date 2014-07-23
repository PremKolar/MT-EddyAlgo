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
    [DD.path.TSow]=tempsalt(DD);
    %% get window according to user input
    [DD.TS.window,~]=GetWindow(DD.path.TSow(1).salt,DD.map.in,DD.TS.keys);
    %% get z info
    DD.TS.window=mergeStruct2(DD.TS.window, GetFields(DD.path.TSow(1).salt, cell2struct({DD.TS.keys.depth},'depth')));
    DD.TS.window.size.Z=numel(DD.TS.window.depth);
    %% distro X lims to chunks
    DD.RossbyStuff.lims.data=limsdata(DD.parameters.RossbySplits,DD.TS.window,'Z');
    %% distro chunks to threads
    DD.RossbyStuff.lims.loop=thread_distro(DD.threads.num,DD.parameters.RossbySplits);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lims=limsdata(splits,window,splitdim)
    %% set dimension for splitting (files dont fit in memory)
    X=window.size.(splitdim);
    %% distro X lims to chunks
    lims=thread_distro(splits,X);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [file]=tempsalt(DD)
    %% find the temp and salt files
    tt=0;ss=0;
    for kk=1:numel(DD.path.TempSalt.files);
        if ~isempty(strfind(upper(DD.path.TempSalt.files(kk).name),DD.TS.keys.salt))
            ss=ss+1;
            file(ss).salt=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name];
        end
        if ~isempty(strfind(upper(DD.path.TempSalt.files(kk).name),DD.TS.keys.temp))
            tt=tt+1;
            file(tt).temp=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name];
        end
    end
end
