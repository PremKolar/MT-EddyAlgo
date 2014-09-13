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
function S01b_fromTS
    %% set up
    %     [DD]=S01b_ST_set_up;
    %   save
    initialise
    load
    %% spmd
%     main(DD)
    %% make netcdf
    WriteNCfile(DD);
    %% update DD
%     save_info(DD);
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
    for cc=lims.loop(id,1):lims.loop(id,2)
        Calculations(DD,cc);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Calculations(DD,cc)
    %% pre-init
    [CK,ccStr]=init(DD,cc,DD.path.Rossby.name);
    if exist(CK.fileSelf,'file') && ~DD.overwrite, return;end
    %% init
    CK=initCK(CK,DD,cc);
    %% calculate Brunt-Väisälä f and potential vorticity
    [CK.N]=calcBrvaPvort(CK,ccStr);
    %% integrate first baroclinic rossby radius
    [CK.rossby.(CK.R1Fname)]=calcRossbyRadius(CK,ccStr);
    %% rossby wave phase speed
    [CK.rossby.(CK.c1Fname)]=calcC_one(CK,ccStr);
    %% clean infs
    CK.N=inf2nan(CK.N);
    for fn=fieldnames(CK.rossby)'
        CK.rossby.(fn{1}) = inf2nan(CK.rossby.(fn{1})) ;
    end
    %% save
    saveChunk(CK);
end
function M=inf2nan(M)
    M(isinf(M))=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CK,ccStr]=init(DD,cc,RossbyDir)
    lims=DD.RossbyStuff.lims.data;
    ccStr=[sprintf(['%0',num2str(length(num2str(size(lims,1)))),'i'],cc),'/',num2str(size(lims,1))];
    disp('initialising..')
    CK.fileSelf=[RossbyDir,'BVRf_',sprintf('%03d',cc),'.mat'];
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
function nc2matSave(DD,fn,in,out,reallocIdx)
    %% get pop data
    in.data=nc_varget(DD.path.Rossby.NCfile,fn);
    data=in.data; %#ok<NASGU>
    %% save pop version either way
    save([DD.path.Rossby.name, fn,'.mat'],'data');clear data;
    if reallocIdx
        %% move pop sized file to another name
        system(['mv ' [DD.path.Rossby.name, fn,'.mat'] ' ' [DD.path.Rossby.name, fn,'PopSize.mat']]);
        data=griddata(in.lon,in.lat,in.data,out.lon,out.lat);
        %% save
        save([DD.path.Rossby.name, fn,'.mat'],'data');
        if strcmp(DD.map.window.type,'globe')
            %% zonal append
            wndw=getfield(load(DD.path.windowFile),'window');
            ovrlpIyx=drop_2d_to_1d(wndw.iy,wndw.ix,size(wndw.iy,1));
            data=data(ovrlpIyx); %#ok<NASGU>
            %% save
            save([DD.path.Rossby.name, fn,'-ZonApp.mat'],'data');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteNCfile(DD)
    splits=DD.parameters.RossbySplits;
    XXlims=DD.RossbyStuff.lims.data;
    yylims=1:DD.TS.window.size.Y;
    
    if any(DD.map.window.fullsize~=DD.TS.window.fullsize)
        reallocIdx=true
    end
    
    %% dummy init
    data=nan([DD.TS.window.size.Y, DD.TS.window.size.X]);  %#ok<NASGU>
    %% loop fields
    for ff=1:numel(DD.FieldKeys.Rossby)
        %% fieldname / fileout name
        FN=DD.FieldKeys.Rossby{ff} ;
        outfileName=[DD.path.Rossby.name FN '.mat'];
        %% start from scratch
        save(outfileName,'data');
        outfile=matfile(outfileName,'Writable',true);
        %% loop chunks
        for cc=1:splits
            xxlims=XXlims(cc,1):XXlims(cc,2);
            CKfn=getfield(getfield(loadChunk(DD.path.Rossby.name,cc),'rossby'),FN);
            outfile.data(yylims,xxlims)=CKfn;
        end
        
        if reallocIdx
            in.lat=nc_varget(DD.path.TempSalt.salt{1},DD.TS.keys.lat);
            in.lon=nc_varget(DD.path.TempSalt.salt{1},DD.TS.keys.lon);
            out.lat=extractdeepfield(load([DD.path.cuts.name DD.path.cuts.files(1).name]),'grids.lat');
            out.lon=extractdeepfield(load([DD.path.cuts.name DD.path.cuts.files(1).name]),'grids.lon');
            
            in.data=load(outfileName,FN)
            %% move pop sized file to another name
            system(['mv ' [outfile] ' ' [DD.path.Rossby.name, FN,'PopSize.mat']]);
           %% resample
            data=griddata(in.lon,in.lat,in.data,out.lon,out.lat);
            %% save
            save([DD.path.Rossby.name, fn,'.mat'],'data');
            
        end
        
      
        
        if strcmp(DD.map.window.type,'globe')
            %% zonal append
            wndw=getfield(load(DD.path.windowFile),'window');
            ovrlpIyx=drop_2d_to_1d(wndw.iy,wndw.ix,size(wndw.iy,1));
            data=data(ovrlpIyx); %#ok<NASGU>
            %% save
            save([DD.path.Rossby.name, fn,'-ZonApp.mat'],'data');
        end
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function nc2mat(DD)
%     %% test for remap
%     reallocIdx=false;
%     fns=ncfieldnames(DD.path.Rossby.NCfile);
%     inTest=numel(nc_varget(DD.path.Rossby.NCfile,fns{1}));
%     if inTest ~= prod(struct2array(DD.map.window.size))
%         reallocIdx=true;
%     end
%     %% save 2 mats
%     [NCin,MATout]=getlatlon(DD);
%     parfor ff=1:numel(fns)
%         fn=fns{ff};
%         nc2matSave(DD,fn,NCin,MATout,reallocIdx);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveChunk(CK)
    save(CK.fileSelf,'-struct','CK');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CK=loadChunk(RossbyDir,cc)
    file_in=[RossbyDir,'BVRf_',sprintf('%03d',cc),'.mat'];
    CK=load(file_in,'rossby');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=	calcRossbyRadius(CK,ccStr)
    dispM(['integrating Rossby Radius for chunk ',ccStr],1)
    [~,YY,XX]=size(CK.N);
    M.depthdiff=repmat(diff(CK.DEPTH),[1 YY XX]);
    % R = 1/(pi f) int N dz
    R=abs(double((squeeze(nansum(M.depthdiff.*CK.N,1))./CK.rossby.f)/pi));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c1]=calcC_one(CK,ccStr)
    %    c=-beta/(k^2+(1/L_r)^2) approx -beta*L^2
    dispM(['applying long rossby wave disp rel for c1 for chunk ',ccStr])
    c1=-CK.rossby.beta.*CK.rossby.(CK.R1Fname).^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N]=calcBrvaPvort(CK,ccStr)
    [ZZ,YY,XX]=size(CK.TEMP);
    dispM(['calculating brunt väisälä, chunk ',ccStr]);
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
function [CK]=initCK(CK,DD,chunk)
    CK.c1Fname=DD.FieldKeys.Rossby{1};
    CK.R1Fname=DD.FieldKeys.Rossby{2};
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
    lat=nc_varget(DD.path.TempSalt.temp{1},DD.TS.keys.lat,dim.start1d, dim.len1d);
    lon=nc_varget(DD.path.TempSalt.temp{1},DD.TS.keys.lon,dim.start1d, dim.len1d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function depth=ChunkDepth(DD)
    depth=nc_varget(DD.path.TempSalt.salt{1},'depth_t');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function salt=ChunkSalt(DD,dim)
    num=numel(DD.path.TempSalt.salt);
    salt=(1/num) * squeeze(nc_varget(DD.path.TempSalt.salt{1},'SALT',dim.start2d,dim.len2d));
    for ss=2:num
        tmp=(1/num) * squeeze(nc_varget(DD.path.TempSalt.salt{ss},'SALT',dim.start2d,dim.len2d));
        salt=salt + tmp;
    end
    salt(salt==0)=nan;
    salt=salt*1000; % to salinity unit. TODO: from input vars
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp=ChunkTemp(DD,dim)
    num=numel(DD.path.TempSalt.temp);
    temp=(1/num) * squeeze(nc_varget(DD.path.TempSalt.temp{1},'TEMP',dim.start2d,dim.len2d));
    for tt=2:num
        tmp=(1/num) * squeeze(nc_varget(DD.path.TempSalt.temp{tt},'TEMP',dim.start2d,dim.len2d));
        temp=temp + tmp;
    end
    temp(temp==0)=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dim=ncArrayDims(DD,cc)
    lims=DD.RossbyStuff.lims.data;
    j_indx_start = DD.TS.window.limits.south-1;
    j_len = DD.TS.window.size.Y;
    dim.start2d = [0 0 j_indx_start lims(cc,1)-1];
    dim.len2d = 	[inf inf j_len diff(lims(cc,:))+1];
    dim.start1d = [j_indx_start lims(cc,1)-1];
    dim.len1d =	[j_len diff(lims(cc,:))+1];
    %% new indeces for output nc file
    xlens=diff(lims,1,2)+1;
    xlens(xlens<0)= xlens(xlens<0) + DD.TS.window.fullsize(2);
    newxstart=sum(xlens(1:cc-1));
    dim.new.start =[0 newxstart];
    dim.new.len =  dim.len1d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


