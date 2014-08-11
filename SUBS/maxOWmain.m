%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 16-Jul-2014 13:52:44
% Computer:GLNX86
% Matlab:7.9
% Author:NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metaData=maxOWmain(DD)
    [Dim,raw] = setup(DD);
    [metaData,sMean]=rhoStuff(DD,raw,Dim);
    calcOW(metaData,raw,sMean);
    save metaD.mat metaData
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [daily,sMean]=rhoStuff(DD,raw,Dim)
    [daily,funcs]=initbuildRho(DD);
    buildRho(daily,raw,Dim,DD.threads.num,funcs) ;
    labBarrier
    %%
    sMean = initbuildRhoMean(DD);
    buildRhoMean(DD.threads.num,sMean,Dim,funcs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcOW(daily,raw,MS)
    T=disp_progress('init','building okubo weiss netcdfs')  ;
    %%
    f=OWfuncs;
    %%
    my=spmdBlockA(MS.Fout,raw,f);
    %%
    toAdd={'OkuboWeiss','log10NegOW'};
    for tt = daily.timesteps;
        T=disp_progress('show',T,numel(daily.timesteps),numel(daily.timesteps));
        loop(daily,f,my,toAdd,tt);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loop(daily,f,my,toAdd,tt)
    OW=fA(f,daily.Fout{tt},daily.OWFout{tt},my,toAdd)      ;
    fB(OW,daily.OWFout{tt},toAdd,f);
end
function OW=fA(f,rhoFile,OWFile,my,toAdd)
    OW=spmdBlockB(f,rhoFile,my);
    initOWNcFile(OWFile,toAdd,size(OW));
end
function fB(OW,OWFile,tA,f)
    f.ncVP(OWFile,OW,tA{1});
    OW(isinf(OW) | OW>=0 | isnan(OW) )=nan;
    f.ncVP(OWFile,log10(-OW),tA{2});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=OWfuncs
    f.ncv=@(d,field) nc_varget(d,field);
    f.ncvOne = @(A) getLocalPart(codistributed(A,codistributor1d(1)));
    f.getHP = @(cf,RM) f.ncvOne(f.ncv(cf,'density')) -  RM;
    f.repinZ = @(A,z) repmat(permute(A,[3,1,2]),[z,1,1]);
    f.ncVP = @(file,OW,field)  nc_varput(file,field,single(OW));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  my=spmdBlockA(MeanFile,raw,f)
    disp('init okubo weiss calcs...')
    spmd
        my.RhoMean=f.ncvOne(f.ncv(MeanFile,'RhoMean'));
        my.Z=size(my.RhoMean,1);
        my.dx=f.repinZ(raw.dx,my.Z);
        my.dy=f.repinZ(raw.dy,my.Z);
        my.GOverF=  f.repinZ(raw.corio.GOverF,my.Z);
        my.depth=f.ncvOne(raw.depth);      
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OW = spmdBlockB(f,currFile,my)
    vc2mstr=@(ow) gcat(ow,1);
    spmd
        my.rhoHighPass=f.getHP(currFile,my.RhoMean);
        ow = vc2mstr(okuweiss(getDefo(UVgrads(getVels(my),my.dx,my.dy))));
        labBarrier
    end
    OW=ow{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ow = okuweiss(d)
    ow = (-(d.vorticity).^2+d.divergence.^2+d.stretch.^2+d.shear.^2)/2;
    %     ow(abs(ow)>1) = nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defo = getDefo(uvg)
    meanin4d = @(A,B) squeeze(mean([permute(A,[4,1,2,3]);permute(B,[4,1,2,3])],1));
    defo.vorticity = uvg.dVdx - uvg.dUdy;
    defo.shear = uvg.dVdx + uvg.dUdy;
    %     defo.divergence = uvg.dUdx + uvg.dVdy;
    %     defo.stretch = uvg.dUdx - uvg.dVdy;
    defo.divergence = 0;
    defo.stretch = - 2* meanin4d(uvg.dVdy,uvg.dUdx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uvg = UVgrads(UV,dx,dy)
    %% calc U gradients
    dUdy = diff(UV.u,1,2);
    dUdx = diff(UV.u,1,3);
    dVdy = diff(UV.v,1,2);
    dVdx = diff(UV.v,1,3);
    uvg.dUdy = dUdy( :,[1:end, end], : )./ dy;
    uvg.dUdx = dUdx( :, :,[1:end, end] )./ dx;
    uvg.dVdy = dVdy( :, [1:end, end], :)./ dy;
    uvg.dVdx = dVdx( :, :,[1:end, end] )./ dx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UV = getVels(my)
    rhoRef = 1000;
    dRho = getDrhodx(my);
    [~,Y,X]=size(my.dx);
    gzOverRhoF = my.GOverF .* repmat(my.depth,[1,Y,X]) / rhoRef;
    UV.u = -dRho.dy .* gzOverRhoF;
    UV.v = dRho.dx .*  gzOverRhoF;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dRho = getDrhodx(my)
    %% calc density gradients
    drdx = diff(my.rhoHighPass,1,3);
    drdy = diff(my.rhoHighPass,1,2);
    dRho.dx = drdx(:,:,[1:end, end]) ./ my.dx;
    dRho.dy = drdy(:,[1:end, end],:) ./ my.dy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = initbuildRhoMean(DD)
    s.files=DD.path.TSow.rho;
    s.Fout=[DD.path.TSow.dailyRhoName 'mean.nc'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,f] = initbuildRho(DD)
    s.timesteps = DD.TSow.lims.timesteps;
    s.keys = DD.TS.keys;
    s.Fin = DD.path.TSow.files;
    s.dirOut=DD.path.full3d.name;
    s.Fout=DD.path.TSow.rho;
    s.OWFout=DD.path.TSow.OW;
    s.geoOut=DD.path.TSow.geo;
    f=buildRhoFuncs;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  buildRhoMean(threads,s,Dim,f)
    initNcFile(s.Fout,'RhoMean',Dim.ws);
    spmd(threads)
        rhoMean = f.locCo(nan(Dim.ws));
        T=disp_progress('init','building density mean')  ;
        labBarrier
    end
    %%
    spmd(threads)
        rhoMean=gop(@vertcat,spmdRhoMeanBlock(f,s,rhoMean,T),1);
        labBarrier
    end
    %%
    f.ncvp(s.Fout,'RhoMean',f.mDit(rhoMean{1},Dim.ws),[0 0 0], [Dim.ws]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhoMean=spmdRhoMeanBlock(f,s,rhoMean,T)
    for ff = 1:numel(s.files)
        T=disp_progress('show',T,numel(s.files),10)  ;
        rhoMean=f.nansumNcvg(rhoMean,s.files{ff},'density');
    end
    rhoMean=rhoMean/numel(s.files);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buildRho(s,raw,Dim,threads,f)
    [depth,lat,T]=spmdInit(threads,raw,Dim,f);
    %%
    for tt = s.timesteps
        if ~exist(s.Fout{tt},'file')
            [RHO,T]=spmdBlock(threads,tt,s,f,T,depth,lat);
            initNcFile(s.Fout{tt},'density',Dim.ws);
            f.ncvp(s.Fout{tt},'density',f.mDit(RHO{1},Dim.ws),[0 0 0], [Dim.ws]);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHO,T]=spmdBlock(threads,tt,s,f,T,depth,lat)
    spmd(threads)
        T=disp_progress('show',T,numel(s.timesteps),numel(s.timesteps))  ;
        RHO=makeRho(s.Fin(tt),depth,lat,s.keys,f)	;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [depth,lat,T]=spmdInit(threads,raw,Dim,f)
    spmd(threads)
        depth = f.locCo(repmat(double(raw.depth),Dim.ws(2)*Dim.ws(3),1));
        lat = f.locCo(f.yx2zyx(raw.lat,Dim.ws(1)));
        T=disp_progress('init','building density netcdfs')  ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=buildRhoFuncs
    f.oneDit = @(md) reshape(md,[],1);
    f.mDit = @(od,ws) reshape(od,[ws(1),ws(2),ws(3)]);
    f.locCo = @(x) getLocalPart(codistributed(f.oneDit(x)));
    f.yx2zyx = @(yx,Z) f.oneDit(repmat(permute(yx,[3,1,2]),[Z,1,1]));
    f.ncvp = @(file,field,array,Ds,De) nc_varput(file,field,array,Ds,De);
    %%
    f.ncvg = @(file,field) nc_varget(file,field);
    f.nansumNcvg = @(A,file,field) nansum([A,f.locCo(f.ncvg(file,field))],2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RHO=makeRho(file,depth,lat,keys,f)
    [temp,salt]=TSget(file,keys,f.locCo);
    pres=f.oneDit(sw_pres(depth,lat));
    RHO=Rhoget(salt,temp,pres);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=Rhoget(salt,temp,pres)
    rho = sw_dens(salt,temp,pres);
    rho(abs(rho>1e10) | rho==0)=nan;
    R=gop(@vertcat, rho,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,S]=TSget(FileIn,keys,locCo)
    T= locCo(nc_varget(FileIn.temp,keys.temp));
    S= locCo(nc_varget(FileIn.salt,keys.salt)*1000);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initOWNcFile(fname,toAdd,WinSize)
    nc_create_empty(fname,'clobber');
    nc_adddim(fname,'k_index',WinSize(1));
    nc_adddim(fname,'i_index',WinSize(3));
    nc_adddim(fname,'j_index',WinSize(2));
    %%
    for kk=1:numel(toAdd)
        ta=toAdd{kk};
        varstruct.Name = ta;
        varstruct.Nctype = 'single';
        varstruct.Dimension = {'k_index','j_index','i_index' };
        nc_addvar(fname,varstruct)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(fname,toAdd,WinSize)
    nc_create_empty(fname,'clobber');
    nc_adddim(fname,'k_index',WinSize(1));
    nc_adddim(fname,'i_index',WinSize(3));
    nc_adddim(fname,'j_index',WinSize(2));
    %%
    varstruct.Name = toAdd;
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'k_index','j_index','i_index' };
    nc_addvar(fname,varstruct)
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

