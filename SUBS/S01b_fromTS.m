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
% NOTE: DOES NOT YET WORK FOR WINDOW SPANNING ZONAL BORDER OF DATA (140611,'yymmdd')

function S01b_fromTS
    %% set up
    [DD]=S01b_ST_set_up;
    %% spmd
    main(DD)
    %% make netcdf
    WriteNCfile(DD);
    %% update DD
    save_info(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD)
    if DD.debugmode
        spmd_body(DD);
    else
        spmd(DD.threads.num)
            spmd_body(DD);
            disp_progress('conclude');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD)
    id=labindex;
    lims=DD.RossbyStuff.lims;
    %% loop over chunks
    for chnk=lims.loop(id,1):lims.loop(id,2)
        Calculations(DD,chnk);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Calculations(DD,chnk)
    %% init
    [CK,cc]=init(DD,chnk);
    %% calculate Brunt-Väisälä f and potential vorticity
    [CK.N]=calcBrvaPvort(CK,cc);
    %% integrate first baroclinic rossby radius
    [CK.rossby.Ro1]=calcRossbyRadius(CK,cc);
    %% rossby wave phase speed
    [CK.rossby.c1]=calcC_one(CK,cc);
    %% save
    disp('saving..')
    saveChunk(CK,DD.path.Rossby.name,chnk);
    %----------------------------------------------------------------------
    function [CK,cc]=init(DD,chnk)
        lims=DD.RossbyStuff.lims.data;
        cc=[sprintf(['%0',num2str(length(num2str(size(lims,1)))),'i'],chnk),'/',num2str(size(lims,1))];
        disp('initialising..')
        CK=initCK(DD,chnk);
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(DD)
    CK=loadChunk(DD.path.Rossby.name,1);
    [dim.Z,dim.Y,dim.X]=size(CK.N);
    nc_adddim(DD.path.Rossby.NCfile,'depth_diff',dim.Z);
    nc_adddim(DD.path.Rossby.NCfile,'i_index',DD.TS.window.size.X);
    nc_adddim(DD.path.Rossby.NCfile,'j_index',DD.TS.window.size.Y);
    %% Ro1
    varstruct.Name = 'RossbyRadius';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'j_index','i_index' };
    nc_addvar(DD.path.Rossby.NCfile,varstruct)
    %% c1
    varstruct.Name = 'RossbyPhaseSpeed';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'j_index','i_index' };
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx,out]=reallocCase(DD)
    %% get lat lon
    [in,out]=getlatlon(DD);
    %% distro all indeces from in to threads
    lims=thread_distro(DD.threads.num,numel(in.lon));
    %% init index vector
    idx=zeros(1,numel(in.lon));
    %% get cross refferencing indeces
    spmd(DD.threads.num)
        JJ=lims(labindex,1):lims(labindex,2);
        idx=getIndicesForOutMaps(in,out,JJ,idx);
        idxx=gop(@vertcat,idx,1);% merge composite
    end
    idx=sum(idxx{1},1);% note idx~=0 only at JJ(labindex). hence vert sum
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [in,out]=getlatlon(DD)
    %% from cuts
    out=load([DD.path.cuts.name DD.path.cuts.files(1).name],'grids');
    out.lat=out.grids.lat;
    out.lon=out.grids.lon;
    out.inc.x=max(mean(diff(out.lon,1,2),1));
    out.inc.y=max(mean(diff(out.lat,1,1),2));
    [out.dim.y,out.dim.x]=size(out.lat);
    out.grids=[];
    %% from pop
    in.lat=nc_varget(DD.path.Rossby.NCfile,'lat');
    in.lon=nc_varget(DD.path.Rossby.NCfile,'lon');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  downscalePop(data,fn,out,DD,idx)
    %% move pop sized file to another name
    system(['mv ' [DD.path.Rossby.name, fn,'.mat'] ' ' [DD.path.Rossby.name, fn,'PopSize.mat']]);
    %% get geo info for output
    [uni,~,~]=unique(idx(idx~=0 & ~isnan(idx)));
    %% nanmean all values out(li) for identical li
    T=disp_progress('init',['remapping ' fn ' from pop data']);
    for li=uni
        T=disp_progress('show',T,numel(uni),10);
        out(li)=nanmean(data(idx==li));
    end
    %% save
    save([DD.path.Rossby.name, fn,'.mat'],'out');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc2matSave(DD,fn,idx,reallocIdx)
    %% get pop data
    in=nc_varget(DD.path.Rossby.NCfile,fn);
    %% save pop version either way
    out=in; %#ok<NASGU>
    save([DD.path.Rossby.name, fn,'.mat'],'out');
    %% remap
    if reallocIdx
        protoout=nan(DD.map.window.size.Y,DD.map.window.size.X);
        downscalePop(in,fn,protoout,DD,idx)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteNCfile(DD)
    %% init
    splits=DD.parameters.RossbySplits;
    initNcFile(DD);
    %% cat chunks
    T=disp_progress('init','creating netcdf');
    for chnk=1:splits
        T=disp_progress('disp',T,splits,100);
        %% put chunks back 2g4
        CK=loadChunk(DD.path.Rossby.name,chnk);
        catChunks2NetCDF(DD.path.Rossby.NCfile,CK);
    end
    %% make also mat files
    nc2mat(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc2mat(DD)
    %% test for remap
    fns=ncfieldnames(DD.path.Rossby.NCfile);
    in=nc_varget(DD.path.Rossby.NCfile,fns{1});
    reallocIdx=false;
    if numel(in)~=prod(struct2array(DD.map.window.size))
        reallocIdx=true;
        [idx,~]=reallocCase(DD);
    else
        idx=[];
    end
    %% save 2 mats
    for ff=1:numel(fns)
        fn=fns{ff};
        nc2matSave(DD,fn,idx,reallocIdx)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function catChunks2NetCDF(file,CK)
    dim=CK.dim.new;
    nc_varput(file,'RossbyRadius',CK.rossby.Ro1,dim.start, dim.len);
    nc_varput(file,'RossbyPhaseSpeed',CK.rossby.c1,dim.start, dim.len);
    nc_varput(file,'lat',CK.lat	,dim.start,dim.len);
    nc_varput(file,'lon',CK.lon	,dim.start, dim.len);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveChunk(CK,RossbyDir,chnk)  %#ok<INUSL>
    file_out=[RossbyDir,'BVRf_',sprintf('%03d',chnk),'.mat'];
    save(file_out,'-struct','CK');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CK=loadChunk(RossbyDir,chnk)
    file_in=[RossbyDir,'BVRf_',sprintf('%03d',chnk),'.mat'];
    CK=load(file_in);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=	calcRossbyRadius(CK,cc)
    dispM(['integrating Rossby Radius for chunk ',cc],1)
    [~,YY,XX]=size(CK.N);
    M.depthdiff=repmat(diff(CK.DEPTH),[1 YY XX]);
    R=abs(double((squeeze(nansum(M.depthdiff.*CK.N,1))./CK.rossby.f)/pi));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c1]=calcC_one(CK,cc)
    %    c=-beta/(k^2+(1/L_r)^2) approx -beta*L^2
    dispM(['applying long rossby wave disp rel for c1 for chunk ',cc])
    c1=-CK.rossby.beta.*CK.rossby.Ro1.^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N]=calcBrvaPvort(CK,cc)
    [ZZ,YY,XX]=size(CK.TEMP);
    dispM(['calculating brunt väisälä, chunk ',cc]);
    %% get full matrices for all variables
    M.depth=double(repmat(CK.DEPTH,[1,YY*XX]));
    M.lat=double(repmat(permute(CK.lat(:),[2 1]), [ZZ,1]));
    M.pressure=double(reshape(sw_pres(M.depth(:),M.lat(:)),[ZZ,YY*XX]));
    M.salt=double(reshape(CK.SALT,[ZZ,YY*XX]));
    M.temp=double(reshape(CK.TEMP,[ZZ,YY*XX]));
    %% get brunt väisälä frequency and pot vort
    [brva,~,~]=sw_bfrq(M.salt,M.temp,M.pressure,M.lat);
    brva(brva<0)=nan;
    N=sqrt(reshape(brva,[ZZ-1,YY,XX]));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CK]=initCK(DD,chunk)
    CK.chunk=chunk;
    CK.dim=ncArrayDims(DD,chunk);
    disp('getting temperature..')
    CK.TEMP=ChunkTemp(DD,CK.dim);
    disp('getting salt..')
    CK.SALT=ChunkSalt(DD,CK.dim);
    disp('getting depth..')
    CK.DEPTH=ChunkDepth(DD);
    disp('getting geo info..')
    [CK.lat,CK.lon]=ChunkLatLon(DD,CK.dim);
    disp('getting coriolis stuff..')
    [CK.rossby]=ChunkRossby(CK);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rossby]=ChunkRossby(CK)
    day_sid=23.9344696*60*60;
    om=2*pi/(day_sid); % frequency earth
    rossby.f=2*om*sind(CK.lat);
    rossby.beta=2*om/earthRadius*cosd(CK.lat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lat,lon]=ChunkLatLon(DD,dim)
    lat=nc_varget(DD.path.TempSalt.temp,DD.TS.keys.lat,dim.start1d, dim.len1d);
    lon=nc_varget(DD.path.TempSalt.temp,DD.TS.keys.lon,dim.start1d, dim.len1d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function depth=ChunkDepth(DD)
	depth=nc_varget(DD.path.TempSalt.salt(1),'depth_t');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function salt=ChunkSalt(DD,dim)
	num=numel(DD.path.TempSalt.salt);
	salt=(1/num) * squeeze(nc_varget(DD.path.TempSalt.salt(1),'TEMP',dim.start2d,dim.len2d));
	for ss=2:num
		tmp=(1/num) * squeeze(nc_varget(DD.path.TempSalt.salt(ss),'TEMP',dim.start2d,dim.len2d));
		salt=salt + tmp;
	end
	salt(salt==0)=nan;
	salt=salt*1000; % to salinity unit. TODO: from input vars
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp=ChunkTemp(DD,dim)
	num=numel(DD.path.TempSalt.temp);
	temp=(1/num) * squeeze(nc_varget(DD.path.TempSalt.temp(1),'TEMP',dim.start2d,dim.len2d));
	for tt=2:num
		tmp=(1/num) * squeeze(nc_varget(DD.path.TempSalt.temp(tt),'TEMP',dim.start2d,dim.len2d));
		temp=temp + tmp;
	end
	temp(temp==0)=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dim=ncArrayDims(DD,chnk)
	lims=DD.RossbyStuff.lims.data;
	j_indx_start = DD.TS.window.limits.south-1;
    j_len = DD.TS.window.size.Y;
    dim.start2d = [0 0 j_indx_start lims(chnk,1)-1];
    dim.len2d = 	[inf inf j_len diff(lims(chnk,:))+1];
    dim.start1d = [j_indx_start lims(chnk,1)-1];
    dim.len1d =	[j_len diff(lims(chnk,:))+1];
    %% new indeces for output nc file
    xlens=diff(squeeze(lims(:,2)));
    xlens(xlens<0)= xlens(xlens<0) + DD.TS.window.fullsize(2);
    newxstart=sum(xlens(1:chnk-1))+1 -1;
    dim.new.start =[0 newxstart];
    dim.new.len =  dim.len1d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


