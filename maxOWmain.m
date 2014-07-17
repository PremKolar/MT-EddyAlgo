%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 16-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOWmain(DD)
    if DD.debugmode
        spmd_body(DD);
    else
        spmd(DD.threads.num)
            spmd_body(DD);
            % 			disp_progress('conclude');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD)
    id=labindex;
    lims=DD.RossbyStuff.lims;
    %%
    [CKpre]=preInitCK(DD);
    %% loop over chunks
    Tf=disp_progress('init','looping through files');
    for ff=0:numel(DD.path.TSow)-1
        Tf=disp_progress('show',Tf,numel(DD.path.TSow),100);
        T=disp_progress('init','looping through files');
        for chnk=lims.loop(id,1):lims.loop(id,2)
            T=disp_progress('yo',T,diff(lims.loop(id))+1,diff(lims.loop(id))+1);
            Calculations(DD,chnk,ff,CKpre);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Calculations(DD,chnk,ff,CK)
    %% init
    lims=DD.RossbyStuff.lims.data;
    cc=[sprintf(['%0',num2str(length(num2str(size(lims,1)))),'i'],chnk),'/',num2str(size(lims,1))];
    dispM('initialising..')
    %% merge
    file_out=[DD.path.Rossby.name,'OW_',sprintf('%03d',ff),'_',sprintf('%03d',chnk),'.mat'];
    if exist(file_out,'file');
        dispM('exists');
        %         return;
    end
    CK=initCK(CK,DD,chnk,ff);
    %% calculate Densisty
    [CK.dens]=calcDens(CK,cc);
    %% OW
    [CK.OW]=calcOW(CK,cc);
    %% save
    dispM('saving..')
    saveChunk(CK,file_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OW =	calcOW(CK,cc)
    dispM(['getting okubo weiss ',cc],1)
    %% rho gradient
    [gr.drdx,gr.drdy]=getDrhodx(CK.dens,CK.DX,CK.DY);
    %% velocities
    vels=getVels(CK.corio,gr,CK.depth);
    clear gr;
    %% uvgrads
    uvg=UVgrads(vels,CK.DX,CK.DY);
    clear vels;
    %% deformation
    deform=getDefo(uvg);
    clear uvg;
    %% okubo weiss
    OW=okuweiss(deform);
    
    
    
    
end
%% ---------------------------------------------------------------------
function OW=okuweiss(d)
    OW=(-d.vorticity.*2+d.divergence.*2+d.stretch.*2+d.shear.*2)/2;
    OW(abs(OW)>1)=nan;
end
%-----------------------------------------------------------------------
function defo=getDefo(uvg)
    defo.vorticity = uvg.dVdx - uvg.dUdy;
    defo.divergence= uvg.dUdx + uvg.dVdy;
    defo.stretch   = uvg.dUdx - uvg.dVdy;
    defo.shear     = uvg.dVdx + uvg.dUdy;
end
%-----------------------------------------------------------------------
function uvg=UVgrads(vels,DX,DY)
    %% calc U gradients
    dUdy=diff(vels.U,1,2);
    dUdx=diff(vels.U,1,3);
    dVdy=diff(vels.V,1,2);
    dVdx=diff(vels.V,1,3);
    [Z,~,~]=size(vels.U);
    uvg.dUdy= dUdy(:, [1:end, end], :)  ./ vertstack(DY,Z);
    uvg.dUdx= dUdx(:, :,[1:end, end] )  ./ vertstack(DX,Z);
    uvg.dVdy= dVdy(:, [1:end, end], :)  ./ vertstack(DY,Z);
    uvg.dVdx= dVdx(:, :,[1:end, end] )  ./ vertstack(DX,Z);
end
%---------
function vels=getVels(cor,gr,depth)
    rhoRef=1000;
    [Z,Y,X]=size(gr.drdy);
    gzOverRhoF=vertstack(cor.GOverF,Z) .* repmat(depth,[1,Y,X]) / rhoRef;
    U=-gr.drdy .* gzOverRhoF;
    V= gr.drdx .* gzOverRhoF;
    
    semi.x=10;
    semi.y=10;
    T=disp_progress('init','high pass filtering geostrophic velocities');
    for z=1:Z
        T=disp_progress('init',T,Z,10);
        [~,vels.U(z,:,:)]=ellipseFltr(semi,squeeze(U(z,:,:)));
        [~,vels.V(z,:,:)]=ellipseFltr(semi,squeeze(V(z,:,:)));
    end
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [drdx,drdy]=getDrhodx(rho,DX,DY)
    %% calc density gradients
    [Z,~,~]=size(rho);
    drdx=diff(rho,1,3);
    drdy=diff(rho,1,2);
    drdx=drdx(:,:,[1:end, end]) ./ vertstack(DX,Z);
    drdy=drdy(:,[1:end, end],:) ./ vertstack(DY,Z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dens]=calcDens(CK,cc)
    [ZZ,YY,XX]=size(CK.TEMP);
    dispM(['calculating pressure, chunk ',cc]);
    %% get full matrices for all variables
    M.depth=double(repmat(CK.depth,[1,YY*XX]));
    M.lat=double(repmat(permute(CK.lat(:),[2 1]), [ZZ,1]));
    
pressure=sw_pres(M.depth(:),M.lat(:)) ; % db 2 Pa
    dens = reshape(sw_dens(CK.SALT(:),CK.TEMP(:),pressure),[ZZ,YY,XX]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CK]=initCK(CK,DD,chunk,ff)
    CK.chunk=chunk;   
    CK.dim=ncArrayDims(DD,chunk,ff);
      CK.depth=ChunkDepth(DD,CK.dim);
    [CK.lat,CK.lon]=ChunkLatLon(DD,CK.dim,ff+1);
    [CK.DY,CK.DX]=ChunkDYDX(CK.lat,CK.lon);
    CK.TEMP=ChunkTemp(DD,CK.dim,ff+1);
    dispM('getting salt..')
    CK.SALT=ChunkSalt(DD,CK.dim,ff+1);
    [CK.rossby]=ChunkRossby(CK);
    CK.corio=coriolisStuff(CK.lat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CK]=preInitCK(DD)
     CK.dim=ncArrayDims(DD,1,0);
     CK.depth=ChunkDepth(DD,CK.dim);
    [CK.lat,CK.lon]=ChunkLatLon(DD,CK.dim,1);
    [CK.DY,CK.DX]=ChunkDYDX(CK.lat,CK.lon);
    [CK.rossby]=ChunkRossby(CK);
    CK.corio=coriolisStuff(CK.lat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=coriolisStuff(lat)
    %% omega
    Omega=angularFreqEarth;
    %% f
    out.f=2*Omega*sind(lat);
    %% beta
    %     out.beta=2*out.Omega/earthRadius*cosd(lat);
    %% gravity
    g=sw_g(lat,zeros(size(lat)));
    %% g/f
    out.GOverF=g./out.f;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rossby]=ChunkRossby(CK)
    day_sid=23.9344696*60*60;
    om=2*pi/(day_sid); % frequency earth
    rossby.f=2*om*sind(CK.lat);
    rossby.beta=2*om/earthRadius*cosd(CK.lat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DY,DX]=ChunkDYDX(lat,lon)
    %% grid increment sizes
    DY=deg2rad(abs(diff(double(lat),1,1)))*earthRadius;
    DX=deg2rad(abs(diff(double(lon),1,2)))*earthRadius.*cosd(lat(:,1:end-1));
    %% append one line/row to have identical size as other fields
    DY=DY([1:end end],:);
    DX=DX(:,[1:end end]);
    %% correct 360° crossings
    seamcrossflag=DX>100*median(DX(:));
    DX(seamcrossflag)=abs(DX(seamcrossflag) - 2*pi*earthRadius.*cosd(lat(seamcrossflag)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveChunk(CK,file_out) %#ok<INUSL>
    save(file_out,'-struct','CK');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dim=ncArrayDims(DD,chnk,ff)
    lims=DD.RossbyStuff.lims.data;
    %% old
    dim.t.old.strt = 0;
    dim.t.old.len  = 1;
    dim.z.old.strt = lims(chnk,1) -1 ;
    dim.z.old.len  = diff(lims(chnk,:)) + 1 ;
    dim.y.old.strt = DD.TS.window.limits.south-1;
    dim.y.old.len  = DD.TS.window.size.Y;
    dim.x.old.strt = DD.TS.window.limits.west-1;
    dim.x.old.len  = DD.TS.window.size.X;   
    %% new indeces for output nc file  
    dim.t.new.strt = ff;
    dim.t.new.len  = 1;    
    dim.z.new.strt = dim.z.old.strt;
    dim.z.new.len  = dim.z.old.len;   
    dim.y.new.strt = 0;
    dim.y.new.len  = DD.TS.window.size.Y;
    dim.x.new.strt = 0;
    dim.x.new.len  = DD.TS.window.size.X;  
    %% array
  [dim]=struct2array(dim);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lat,lon]=ChunkLatLon(DD,dim,ff)
   
    lat=nc_varget(DD.path.TSow(ff).temp,DD.TS.keys.lat...
        ,edf(dim,'old.strt',3,4),edf(dim,'old.len',3,4));
    
    lon=nc_varget(DD.path.TSow(ff).temp,DD.TS.keys.lon...
        ,edf(dim,'old.strt',3,4),edf(dim,'old.len',3,4));       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function depth=ChunkDepth(DD,dim)
    depth=nc_varget(DD.path.TSow(1).salt,DD.TS.keys.depth...
       ,edf(dim,'old.strt',2,2),edf(dim,'old.len',2,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function salt=ChunkSalt(DD,dim,ff)
    salt=squeeze(nc_varget(DD.path.TSow(ff).salt,DD.TS.keys.salt...
        ,edf(dim,'old.strt'),edf(dim,'old.len')));
  
    salt(salt==0)=nan;
    salt=salt*1000; % to salinity unit. TODO: from input vars
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp=ChunkTemp(DD,dim,ff)
    temp=squeeze(nc_varget(DD.path.TSow(ff).temp,DD.TS.keys.temp...
       ,edf(dim,'old.strt'),edf(dim,'old.len')));
  
   temp(temp==0)=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% shorthand
function out=edf(in,field,from,till)
    out=extractdeepfield(in,field);
    if nargin == 4
        out=out(from:till);
    end
end
