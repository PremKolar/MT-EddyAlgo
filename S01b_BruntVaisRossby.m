%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% needs one 3D salt and temperature file each.
% integrates over depth to calculate
% -Brunt Väisälä frequency
% -Rossby Radius
% -Rossby wave first baroclinic phase speed
% -Pot. vorticity
function S01b_BruntVaisRossby
    %% set up
    [DD,lims]=set_up;
    %% spmd
    spmd(DD.threads.num)
        spmd_body(DD,lims);
    end
    %% make netcdf
    DD.path.Rossby=WriteNCfile(DD,lims);
    %% update DD
    save_info(DD);
end
function [DD,lims]=set_up
    %% init
    DD=initialise;
    %% set number of chunks to split large data (see input_vars.m)
    splits=DD.RossbyStuff.splits;
    %% threads
    DD.threads.num=init_threads(DD.threads.num);
    %% find temp and salt files
    [DD.path.TempSalt.salt,DD.path.TempSalt.temp]=tempsalt(DD);
    %% set dimension for splitting (files dont fit in memory)
    X=DD.map.window.size.X;
    %% map chunks
    lims.data=thread_distro(splits,X) + DD.map.window.limits.west-1;
    %% distro chunks to threads
    lims.loop=thread_distro(DD.threads.num,splits);
end
function spmd_body(DD,lims)
    id=labindex;
    %% loop over chunks
    for chnk=lims.loop(id,1):lims.loop(id,2)
        Calculations(DD,lims.data,chnk);
    end
end
function Calculations(DD,lims,chnk)
    cc=[sprintf(['%0',num2str(length(num2str(size(lims,1)))),'i'],chnk),'/',num2str(size(lims,1))];
    disp('initialising..')
    CK=initCK(DD,lims,chnk);
    %% calculate Brunt-Väisälä f and potential vorticity
    [CK.BRVA,CK.PVORT,CK.midDepth]=calcBrvaPvort(CK,cc);
    %% integrate first baroclinic rossby wave phase speed
    [CK.rossby.c1]=calcC_one(CK,cc);
    %% rossby radius ie c1/f
    [CK.rossby.Ro1]=calcRossbyRadius(CK.rossby);
    %% save
    disp('saving..')
    saveChunk(CK,DD,chnk);
end
function nc_file_name = WriteNCfile(DD,lims)
    nc_file_name=initNC(DD);
    splits=DD.RossbyStuff.splits;
    T=disp_progress('init','creating netcdf');
    for chnk=1:splits
        T=disp_progress('disp',T,splits,100);
        %% put chunks back 2g4
        catChunks2NetCDF(DD,lims,chnk,nc_file_name)
    end
end
function	nc_file_name = initNC(DD)
    nc_file_name=[DD.path.Rossby.name, 'BVRf_all.nc'];
    overwriteornot(nc_file_name)
end

function overwriteornot(nc_file_name)
    %% test if file exists already
    try
        nc_create_empty(nc_file_name,'noclobber')
    catch me
        disp(me.message)
        reply = input('Do you want to overwrite? Y/N [Y]: ', 's');
        if isempty(reply)
            reply = 'Y';
        end
        if strcmp(reply,'Y')
            nc_create_empty(nc_file_name,'clobber')
        else
            error('exiting')
        end
    end
end

function catChunks2NetCDF(DD,lims,chnk,nc_file_name)
    %% init
    CK=loadChunk(DD,chnk);
    [Z,Y,X]=size(CK.rossby.c1);
    strt=lims.data(chnk,1)-DD.map.window.limits.west;
    dim.fourD.start  = [ 0 0 strt];
    dim.fourD.length = [Z Y X];
    dim.twoD.start   = [0];
    dim.twoD.length = [Z];
    %% add dimensions to netcdf
    if chnk==1
        nc_adddim(nc_file_name,'depth_diff',Z);
        nc_adddim(nc_file_name,'i_index',DD.map.window.size.X	);
        nc_adddim(nc_file_name,'j_index',DD.map.window.size.Y);
    end
    %% N
    varstruct.Name = 'BruntVaisala';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'depth_diff','j_index','i_index' };
    if chnk==1, nc_addvar(nc_file_name,varstruct); end
    nc_varput(nc_file_name,'BruntVaisala',CK.BRVA,dim.fourD.start, dim.fourD.length);
    %% Ro1
    varstruct.Name = 'RossbyRadius';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'depth_diff','j_index','i_index' };
    if chnk==1,nc_addvar(nc_file_name,varstruct); end
    nc_varput(nc_file_name,'RossbyRadius',CK.rossby.Ro1,dim.fourD.start, dim.fourD.length);
    %% c1
    varstruct.Name = 'RossbyPhaseSpeed';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'depth_diff','j_index','i_index' };
    if chnk==1,nc_addvar(nc_file_name,varstruct); end
    nc_varput(nc_file_name,'RossbyPhaseSpeed',CK.rossby.c1,dim.fourD.start, dim.fourD.length);
    %% depth
    if chnk==1
        varstruct.Name = 'Depth';
        varstruct.Nctype = 'double';
        varstruct.Dimension ={'depth_diff'};
        nc_addvar(nc_file_name,varstruct)
        nc_varput(nc_file_name,'Depth',CK.midDepth,dim.twoD.start, dim.twoD.length);
    end
end
function saveChunk(CK,DD,chnk) %#ok<*INUSL>
    file_out=[DD.path.TempSalt.name,'BVRf_',sprintf('%03d',chnk),'.mat'];
    save(file_out,'-struct','CK');
end
function CK=loadChunk(DD,chnk)
    file_in=[DD.path.TempSalt.name,'BVRf_',sprintf('%03d',chnk),'.mat'];
    CK=load(file_in);
end
function R=	calcRossbyRadius(rossby)
    %% lambda=c/f
    f=repmat(permute(rossby.f,[3 1 2]),[size(rossby.c1,1) 1 1]);
    R=abs(rossby.c1./f);
end
function [c1]=calcC_one(CK,cc)
    c1=CK.rossby.c1;
    [ZZ,YY,XX]=size(c1);
    %% int_-h^0 N/pi dz
    disp(['integrating phase speed for chunk ',cc])
    M.depthdiff=repmat(diff(CK.DEPTH),[1 YY XX]);
    for zz=1:ZZ
        c1(zz,:,:)=-squeeze(nansum(M.depthdiff(1:zz,:,:).*CK.BRVA(1:zz,:,:),1))/pi;
    end
    
end
function [BRVA,PVORT,midDepth]=calcBrvaPvort(CK,cc)
    [ZZ,YY,XX]=size(CK.BRVA);
    disp(['calculating brunt väisälä, chunk ',cc]);
    %% get full matrices for all variables
    M.depth=repmat(CK.DEPTH,[1,YY*XX]);
    M.lat=repmat(permute(CK.LAT(:),[2 1]), [ZZ+1,1]);
    M.pressure=reshape(sw_pres(M.depth(:),M.lat(:)),[ZZ+1,YY*XX]);
    M.salt=reshape(CK.SALT,[ZZ+1,YY*XX]);
    M.temp=reshape(CK.TEMP,[ZZ+1,YY*XX]);
    %% get brunt väisälä frequency and pot vort
    [brva,pvort,md]=sw_bfrq(M.salt,M.temp,M.pressure,M.lat);
    midDepth=squeeze(md(:,1)); % uniform
    BRVA=sqrt(reshape(brva,[ZZ,YY,XX]));
    BRVA(abs(imag(BRVA))>0)=0;
    PVORT=reshape(pvort,[ZZ,YY,XX]);
    
end
function [CK,DD]=initCK(DD,lims,chnk)
    CK.chunk=chnk;
    dim=ncArrayDims(DD,lims,chnk);
    disp('getting temperature..')
    CK.TEMP=ChunkTemp(DD,dim);
    disp('getting salt..')
    CK.SALT=ChunkSalt(DD,dim);
    disp('getting depth..')
    CK.DEPTH=ChunkDepth(DD);
    disp('getting geo info..')
    [CK.LAT,CK.LON]=ChunkLatLon(DD,dim);
    disp('init pot vort and N..')
    [CK.BRVA,CK.PVORT]=ChunkBrvaPvort(size(CK.TEMP));
    disp('init rossby radius/phase speed..')
    [CK.rossby]=ChunkRossby(CK);
end
function [rossby]=ChunkRossby(CK)
    day_sid=23.9344696*60*60;
    om=2*pi/(day_sid); % frequency earth
    rossby.Ro1=nan(size(CK.BRVA));
    rossby.c1=nan(size(CK.BRVA));
    rossby.f=2*om*sind(CK.LAT);
end
function [BRVA,PVORT]=ChunkBrvaPvort(sze)
    BRVA=nan(sze);
    PVORT=nan(sze);
    %% append one more level so data is of equal dims
    BRVA(end,:,:)=[];
    PVORT(end,:,:)=[];
end
function [lat,lon]=ChunkLatLon(DD,dim)
    % TODO: get nc var strings from input_vars.m
    lat=nc_varget(DD.path.TempSalt.temp,'U_LAT_2D',dim.twoD.start, dim.twoD.length);
    lon=nc_varget(DD.path.TempSalt.temp,'U_LON_2D',dim.twoD.start, dim.twoD.length);
end
function depth=ChunkDepth(DD)
    depth=nc_varget(DD.path.TempSalt.salt,'depth_t');
end
function salt=ChunkSalt(DD,dim)
    dispNcInfo(DD.path.TempSalt.salt)
    salt=squeeze(nc_varget(DD.path.TempSalt.salt,'SALT',dim.fourD.start,dim.fourD.length));
    salt(salt==0)=nan;
    salt=salt*1000; % to salinity unit. TODO: from input vars
end
function temp=ChunkTemp(DD,dim)
    dispNcInfo(DD.path.TempSalt.temp)
    temp=squeeze(nc_varget(DD.path.TempSalt.temp,'TEMP',dim.fourD.start,dim.fourD.length));
    temp(temp==0)=nan;
end
function dispNcInfo(ncIn)
    %% works for the POP data...
    try
        info=nc_info(ncIn);
        disp(info.Dataset(end-1).Name);
        disp(info.Dataset(end-1).Dimension);
    catch disperr
        disp(disperr)
    end
end
function dim=ncArrayDims(DD,lims,chnk)
    j_indx_start = DD.map.window.limits.south;
    j_len = DD.map.window.size.Y;
    dim.fourD.start = [0 0 j_indx_start lims(chnk,1)-1];
    dim.fourD.length = 	[inf inf j_len diff(lims(chnk,:))+1];
    dim.twoD.start = [j_indx_start lims(chnk,1)-1];
    dim.twoD.length =	[j_len diff(lims(chnk,:))+1];
end
function [fileS,fileT]=tempsalt(DD)
    %% find the temp and salt files
    for kk=1:numel(DD.path.TempSalt.files)
        if ~isempty(strfind(DD.path.TempSalt.files(kk).name,'SALT'))
            fileS=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name];
        end
        if ~isempty(strfind(DD.path.TempSalt.files(kk).name,'TEMP'))
            fileT=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name];
        end
    end
end

