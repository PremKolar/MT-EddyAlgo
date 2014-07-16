%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 16-Jul-2014 13:53:44
% Computer:GLNX86
% Matlab:7.9
% Author:NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOWwriteNCfile(DD)
    %% init
    splits=DD.parameters.RossbySplits;
    initNcFile(DD);
    Tf=disp_progress('init','looping over time-steps');
    fnum=numel(DD.path.TSow);
    for ff=0:fnum-1
        Tf=disp_progress('disp',Tf,fnum,100);
        %% cat chunks
        T=disp_progress('init','creating netcdf');
        for chnk=1:splits
            T=disp_progress('disp',T,splits,6);
            %% put chunks back 2g4
            CK=loadChunk(DD.path.Rossby.name,chnk,ff);
            catChunks2NetCDF(DD.path.Rossby.NCfile,CK,ff);
        end
    end
    %% make also mat files
    nc2mat(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CK=loadChunk(RossbyDir,chnk,ff)
    file_in=[RossbyDir,'OW_',sprintf('%03d',ff),'_',sprintf('%03d',chnk),'.mat'];
    CK=load(file_in);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(DD)
    dim.t=numel(DD.path.TSow);
    nc_adddim(DD.path.Rossby.NCfile,'k_index',DD.TS.window.size.Z);
    nc_adddim(DD.path.Rossby.NCfile,'i_index',DD.TS.window.size.X);
    nc_adddim(DD.path.Rossby.NCfile,'j_index',DD.TS.window.size.Y);
    nc_adddim(DD.path.Rossby.NCfile,'t_index',dim.t);
    %% OW
    varstruct.Name = 'OkuboWeiss';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'t_index','k_index','j_index','i_index' };
    nc_addvar(DD.path.Rossby.NCfile,varstruct)
    %% depth
    varstruct.Name = 'depth';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'k_index'};
    nc_addvar(DD.path.Rossby.NCfile,varstruct)
    %% lat
    varstruct.Name = 'lat';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'j_index','i_index' };
    nc_addvar(DD.path.Rossby.NCfile,varstruct)
    %% lon
    varstruct.Name = 'lon';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'j_index','i_index' };
    nc_addvar(DD.path.Rossby.NCfile,varstruct)
    %% time
    varstruct.Name = 'time';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'t_index'};
    nc_addvar(DD.path.Rossby.NCfile,varstruct)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function catChunks2NetCDF(file,CK,ff)
    start=extractdeepfield(CK,'dim.new.strt');
    len=extractdeepfield(CK,'dim.new.len');
    
    nc_varput(file,'OkuboWeiss',CK.OW, start,len);
    nc_varput(file,'depth',CK.depth,start(2),len(2));
    nc_varput(file,'lat', CK.lat,start(3:4),len(3:4));
    nc_varput(file,'lon', CK.lon,start(3:4),len(3:4));
    nc_varput(file,'time', now + ff,[ff], [1]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc2mat(DD)
    %% test for remap
    fns=ncfieldnames(DD.path.Rossby.NCfile)';
    %% save 2 mats
    files=struct;
    for fn=fns;
        out=nc_varget(DD.path.Rossby.NCfile,fn{1});
        files.(fn{1})=nc2matSave(DD,out,fn{1});
    end
    %% save mat file info
    save([DD.path.Rossby.name 'filenames.mat'],'files')
   
end

 %----------------------------------------------------------------------
    function filename=nc2matSave(DD,out,fn) %#ok<INUSL>
        filename=[DD.path.Rossby.name, fn,'.mat'];
        save([DD.path.Rossby.name, fn,'.mat'],'out');
    end

