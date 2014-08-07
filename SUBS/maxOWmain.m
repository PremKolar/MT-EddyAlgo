%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 16-Jul-2014 13:52:44
% Computer:GLNX86
% Matlab:7.9
% Author:NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metaD=maxOWmain(DD)
%     if exist('metaD.mat','file')
%         load metaD
%         return
%     end
    [Dim,raw] = setup(DD);
    if DD.debugmode
        d=spmd_body(DD,Dim,raw);
    else
        spmd(DD.threads.num)
            d=spmd_body(DD,Dim,raw);
        end
    end
    metaD=d{1};
    metaD.dim=Dim;
    
    save metaD.mat metaD
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=spmd_body(DD,Dim,raw)
     d.daily=initbuildRho(Dim,DD);
     buildRho(d.daily,raw,Dim) ;
    labBarrier
    [s] = initbuildRhoMean(DD);
    buildRhoMean(s,Dim)
    labBarrier
    calcOW(d,Dim,raw,s);
    labBarrier
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcOW(d,Dim,raw,MeanStuff)
    
    todo=d.daily.timesteps;
    T=disp_progress('init','building okubo weiss netcdfs')  ;
    for zz = 0:Dim.ws(1)-1
        T=disp_progress('show',T,Dim.ws(1),Dim.ws(1))  ;
        strt = [zz 0 0] ;
        len = [1 Dim.ws(2:3)] ;
        Tt=disp_progress('init','looping over time')  ;
        for tt = todo;
            Tt=disp_progress('show',Tt,numel(d.daily.timesteps),10)  ;
            rhoHighPass=nc_varget(d.daily.Fout{tt+1},'density',strt,len) - nc_varget(MeanStuff.Fout,'RhoMean',strt,len) ;
            %%
            %             depth = nc_varget(d.daily.geoOut,'depth');
            %% rho gradient
            [gr.drdx,gr.drdy] = getDrhodx(rhoHighPass,raw.dx,raw.dy);
            %% velocities
            vels = getVels(raw.corio.GOverF,gr,raw.depth(zz+1));
            %% uvgrads
            uvg = UVgrads(vels,raw.dx,raw.dy);
            %% deformation
            deform = getDefo(uvg);
            %% okubo weiss
            OW = permute(okuweiss(deform),[3,1,2]);
            %             nc_varput(d.daily.OWFout{tt+1},'OkuboWeiss',OW,strt,len);
            NCcritWrite(d.daily.OWFout{tt+1},'OkuboWeiss',OW,strt,len)
        end
    end
    
    
    
    %%
    %     for tt = todo;
    %         system(['mv ' d.daily.OWFout{tt+1} ' '  finishedF{tt+1} ]);
    %     end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ow = okuweiss(d)
    ow = (-d.vorticity.*2+d.divergence.*2+d.stretch.*2+d.shear.*2)/2;
    ow(abs(ow)>1) = nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defo = getDefo(uvg)
    defo.vorticity = uvg.dVdx - uvg.dUdy;
    defo.divergence = uvg.dUdx + uvg.dVdy;
    defo.stretch = uvg.dUdx - uvg.dVdy;
    defo.shear = uvg.dVdx + uvg.dUdy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uvg = UVgrads(vels,dx,dy)
    %% calc U gradients
    dUdy = diff(vels.U,1,1);
    dUdx = diff(vels.U,1,2);
    dVdy = diff(vels.V,1,1);
    dVdx = diff(vels.V,1,2);
    uvg.dUdy = dUdy( [1:end, end], :)./ dy;
    uvg.dUdx = dUdx( :,[1:end, end] )./ dx;
    uvg.dVdy = dVdy( [1:end, end], :)./ dy;
    uvg.dVdx = dVdx( :,[1:end, end] )./ dx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vels = getVels(GOverF,gr,depth)
    rhoRef = 1000;
    gzOverRhoF = GOverF * depth / rhoRef;
    vels.U = -gr.drdy .* gzOverRhoF;
    vels.V = gr.drdx .* gzOverRhoF;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [drdx,drdy] = getDrhodx(rho,dx,dy)
    %% calc density gradients
    drdx = diff(rho,1,2);
    drdy = diff(rho,1,1);
    drdx = drdx(:,[1:end, end]) ./ dx;
    drdy = drdy([1:end, end],:) ./ dy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = initbuildRhoMean(DD)
    
    s.lims = DD.TSow.lims.inZ-1;
    s.Zsteps = s.lims(labindex,1):s.lims(labindex,2);
    s.files=DD.path.TSow.rho;
    s.Fout=[DD.path.TSow.dailyRhoName 'mean.nc'];
    initNcFile(s.Fout,'RhoMean',DD.TSow.window.size);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = initbuildRho(Dim,DD)
    dispM('building Rho',1)
    s.id = labindex;
    s.lims = DD.TSow.lims.inTime-1;
    s.timesteps = s.lims(s.id,1):s.lims(s.id,2);
    s.keys = DD.TS.keys;
    s.Ysteps=round(linspace(0,Dim.ws(2)-1,100));
    T=disp_progress('init','init Rho')  ;
    for cc =1:numel(s.Ysteps)-1
        T=disp_progress('show',T,numel(s.Ysteps)-1,100)  ;
        s.Ychunks{cc}=s.Ysteps(cc):s.Ysteps(cc+1)-1;
    end
    s.Fin = DD.path.TSow.files;
    s.dirOut=DD.path.full3d.name;
    s.Fout=DD.path.TSow.rho;
    s.OWFout=DD.path.TSow.OW;
    s.geoOut=DD.path.TSow.geo;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  buildRhoMean(s,Dim)
    T=disp_progress('init','building density mean')  ;
    myzz=0;
    for zz = s.Zsteps
        jhvuvgh
        myzz=myzz+1;
        T=disp_progress('show',T,numel(s.Zsteps),10)  ;
        strt=[zz 0 0 ];
        len =[1  Dim.ws(2:3)  ];
        rhoAtZ=meanOverTime(s,Dim,strt,len);
        rhoAtZ=reshape(rhoAtZ/numel(s.files),Dim.ws(2),[]);
        rhoAtZ(rhoAtZ>1e10)=nan;
        NCcritWrite(s.Fout,'RhoMean',rhoAtZ,strt,len)
    end
end
%----------------------------------------------------------------------
function rhoAtZ=meanOverTime(s,Dim,strt,len)
    rhoAtZ=nan(Dim.ws(2)*Dim.ws(3),1);
    Tr=disp_progress('init','looping over files')  ;
    for ff = 1:numel(s.files)
        Tr=disp_progress('show',Tr,numel(s.files),10)  ;
        rhoAtZ=nansum([rhoAtZ, reshape(nc_varget(s.files{ff},'density',strt,len),[],1)],2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buildRho(s,raw,Dim)
    Dim.strt=[0 0 0];
    Dim.len=Dim.ws;
    DimOri=Dim;
    T=disp_progress('init','building density netcdfs')  ;
    for tt = s.timesteps
        T=disp_progress('show',T,numel(s.timesteps),10)  ;
        Dim=DimOri;
        Ty=disp_progress('init','looping over y chunks')  ;
        initNcFile(s.Fout{tt+1},'density',Dim.ws)
        for cc = 1:numel(s.Ysteps)-1
            Ty=disp_progress('show',Ty,numel(s.Ysteps),5)  ;
            Dim.strt(2)  = DimOri.strt(2) + s.Ysteps(cc);
            Dim.len(2)  = length(s.Ychunks{cc});
            [TS] = getNowAtY(s.Fin(tt+1),s.keys,Dim,s.dirOut,raw,s.Ychunks{cc}+1);
            rho = reshape(calcDens(TS,Dim.len),Dim.len);
            nc_varput(s.Fout{tt+1},'density',rho,Dim.strt,Dim.len);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(fname,toAdd,WinSize)
    nc_create_empty(fname,'clobber');
    nc_adddim(fname,'k_index',WinSize(1));
    nc_adddim(fname,'i_index',WinSize(3));
    nc_adddim(fname,'j_index',WinSize(2));
    %% rho
    varstruct.Name = toAdd;
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'k_index','j_index','i_index' };
    nc_addvar(fname,varstruct)
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = getNowAtY(FileIn,keys,Dim,dirOut,raw,YY)
    make1d = @(Axb,C) reshape(repmat(double(permute(Axb,[3,1,2])),C,1),[],1)   ;
    dsr = @(M) double(sparse(reshape(M,[],1)));
    Z=Dim.ws(1);
    out.lat = make1d(raw.lat(YY,:),Z);
    out.lon = make1d(raw.lon(YY,:),Z);
    out.dx = make1d(raw.dx(YY,:),Z);
    out.dy = make1d(raw.dy(YY,:),Z);
    out.depth = repmat(double(raw.depth),numel(YY),1);
    out.rawDepth = raw.depth;
    out.temp = dsr(nc_varget(FileIn.temp,keys.temp,[0 Dim.strt],[1 Dim.len]));
    out.salt = dsr(nc_varget(FileIn.salt,keys.salt,[0 Dim.strt],[1 Dim.len])*1000);
    out.file = fileStuff(FileIn.salt, dirOut);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHO] = calcDens(TS,Dim)
    depth=repmat(TS.depth,Dim(3),1);
    RHO=dens(TS,depth)   ;
    function rho=dens(TS,depth)
        rho=sw_dens(TS.salt,TS.temp,sw_pres(depth,TS.lat));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dim,raw] = setup(DD)
    ws = DD.TSow.window.size;
    FileIn = DD.path.TSow.files(1);
    keys = DD.TS.keys;
    raw.depth = nc_varget(FileIn.salt,keys.depth);
    raw.lat = nc_varget(FileIn.temp, keys.lat);
    raw.lon = nc_varget(FileIn.temp, keys.lon);
    [raw.dy,raw.dx] = getdydx( raw.lat, raw.lon);
    raw.corio = coriolisStuff(raw.lat);
    Dim.ws=ws;
    
    %% geo
    nc_varput(DD.path.TSow.geo,'depth',raw.depth);
    nc_varput(DD.path.TSow.geo,'lat',raw.lat);
    nc_varput(DD.path.TSow.geo,'lon',raw.lon);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = fileStuff(fin,dirout)
    dateIdx = regexp(fin,'[0-9]{8}');
    dateN = fin(dateIdx:dateIdx+7);
    out.date = datenum(dateN,'yyyymmdd');
    out.fout = [ dirout, 'meanDens',dateN,'.mat'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = coriolisStuff(lat)
    OmegaTw = 2*angularFreqEarth;
    %% f
    out.f = OmegaTw*sind(lat);
    %% beta
    out.beta = OmegaTw/earthRadius*cosd(lat);
    %% gravity
    g = sw_g(lat,zeros(size(lat)));
    %% g/f
    out.GOverF = g./out.f;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dy,dx] = getdydx(lat,lon)
    %% grid increment sizes
    dy = deg2rad(abs(diff(double(lat),1,1)))*earthRadius;
    dx = deg2rad(abs(diff(double(lon),1,2)))*earthRadius.*cosd(lat(:,1:end-1));
    %% append one line/row to have identical size as other fields
    dy = dy([1:end end],:);
    dx = dx(:,[1:end end]);
    %% correct 360° crossings
    seamcrossflag = dx>100*median(dx(:));
    dx(seamcrossflag) = abs(dx(seamcrossflag) - 2*pi*earthRadius.*cosd(lat(seamcrossflag)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

