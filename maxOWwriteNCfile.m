
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 16-Jul-2014 13:53:44
% Computer:GLNX86
% Matlab:7.9
% Author:NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DD=maxOWwriteNCfile(DD)
       [dirstuff]=init(DD);
    %% daily operations
    for ff=0:dirstuff.fnum-1
        %% init        
       dayloop(ff,dirstuff,DD.TS.window.size);
    end
    %% make also large mat file
%     nc2mat(DD);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tf]=dayloop(ff,dirstuff,WinSize)
    Tf=disp_progress('disp',Tf,dirstuff.fnum,100);
    dirstuff.fname=initNcFile(WinSize,ff,dirstuff);   
    NCLoop(dirstuff,ff);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NCLoop(dirstuff,ff)
    dispM('creating netcdf');
    do=dirstuff.out;
    dfn=dirstuff.fname;
    parfor chnk=1:dirstuff.splits
        %% put chunks back 2g4
        backpackNC(do,dfn,ff,chnk);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function backpackNC(do,dfn,ff,chnk)
    CK=loadChunk(do,chnk,ff);
    catChunks2NetCDF(dfn,CK,ff);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dirstuff]=init(DD)    
    dirstuff.fnum=numel(DD.path.TSow);
    dirstuff.out   = [DD.path.Rossby.name];  
dirstuff.dailyBaseName   = [ dirstuff.out 'OW_']; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CK=loadChunk(RossbyDir,chnk,ff)
    file_in=[RossbyDir,'OW_',sprintf('%03d',ff),'_',sprintf('%03d',chnk),'.mat'];
    CK=load(file_in);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fname=initNcFile(WinSize,ff,ds)
   fname=[ds.dailyBaseName sprintf('%04d.nc',ff) ];
   nc_create_empty(fname,'clobber');
    nc_adddim(fname,'k_index',WinSize.Z);
    nc_adddim(fname,'i_index',WinSize.X);
    nc_adddim(fname,'j_index',WinSize.Y);
%     nc_adddim(fname,'t_index',dim.t);
    nc_adddim(fname,'t_index',1);
    %% OW
    varstruct.Name = 'OkuboWeiss';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'t_index','k_index','j_index','i_index' };
    nc_addvar(fname,varstruct)
    %% depth
    varstruct.Name = 'depth';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'k_index'};
    nc_addvar(fname,varstruct)
    %% lat
    varstruct.Name = 'lat';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'j_index','i_index' };
    nc_addvar(fname,varstruct)
    %% lon
    varstruct.Name = 'lon';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'j_index','i_index' };
    nc_addvar(fname,varstruct)
    %% time
    varstruct.Name = 'time';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'t_index'};
    nc_addvar(fname,varstruct)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function catChunks2NetCDF(file,CK,timestep)
    start=extractdeepfield(CK,'dim.new.strt');
     start(1) = 0;
    len=extractdeepfield(CK,'dim.new.len');
    
    nc_varput(file,'OkuboWeiss',CK.OW, start,len);
    nc_varput(file,'depth',CK.depth,start(2),len(2));
    nc_varput(file,'lat', CK.lat,start(3:4),len(3:4));
    nc_varput(file,'lon', CK.lon,start(3:4),len(3:4));
%     nc_varput(file,'time', now + timestep,[timestep], [1]);
    nc_varput(file,'time', timestep       ,[0], [1]);
end

% function nc2mat(DD)
%     %% test for remap
%     fns=ncfieldnames(DD.path.Rossby.NCfile)';
%     %% save 2 mats
%     files=struct;
%     for fn=fns;
%         out=nc_varget(DD.path.Rossby.NCfile,fn{1});
%         files.(fn{1})=nc2matSave(DD,out,fn{1});
%     end
%     %% save mat file info
%     save([DD.path.Rossby.name 'filenames.mat'],'files')
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function filename=nc2matSave(DD,out,fn) %#ok<INUSL>
%     filename=[DD.path.Rossby.name, fn,'.mat'];
%     if exist(filename,'file'), return; end
%     save([DD.path.Rossby.name, fn,'.mat'],'out');
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function saveMatChunks(outDir,dim,OW,timeStep)
%     [Z,Y,X]=size(OW);
%     start = extractdeepfield(dim,'new.strt') + 1;
%     len   = extractdeepfield(dim,'new.len') ;
%     till  = start + len -1 ;
%     chunkFname = [ outDir 'OWchunk' sprintf('%04dtime_%04d-%04ddepth.mat',timeStep,start(2),till(2)) ];
%     OW2d=reshape(permute(OW,[2,3,1]),[Y,X*Z]); %#ok<NASGU>
%     save(chunkFname,'OW2d');   
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%