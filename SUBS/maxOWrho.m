%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 16-Jul-2014 13:52:44
% Computer:GLNX86
% Matlab:7.9
% Author:NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOWrho
    load DD
    DD.MD=main(DD,DD.raw,DD.Dim,DD.f);
    save DD
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MD]=main(DD,raw,Dim,f)
     [MD]=initbuildRho(DD);
    buildRho(MD,raw,Dim,DD.threads.num,f) ;
    MD.sMean = initbuildRhoMean(DD);
    buildRhoMean(DD.threads.num,MD.sMean,Dim,f);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = initbuildRhoMean(DD)
    s.files=DD.path.TSow.rho;
    s.Fout=[DD.path.TSow.dailyRhoName 'mean.nc'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = initbuildRho(DD)
   dF
   s.timesteps= DD.TSow.lims.timesteps;
    s.keys = DD.TS.keys;
    s.Fin = DD.path.TSow.files;
    s.dirOut=DD.path.full3d.name;
    s.Fout=DD.path.TSow.rho;
    s.OWFout=DD.path.TSow.OW;
    s.geoOut=DD.path.TSow.geo;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  buildRhoMean(threads,s,Dim,f)
    if ~exist(s.Fout,'file')
        selMstr=@(x) x{1};
        rhoMean=selMstr(buildRhoMeanOperate(threads,s,Dim,f));
        f.ncvp([s.Fout 'tmp'],'RhoMean',f.mDit(rhoMean,Dim.ws),[0 0 0], [Dim.ws]);
        system(['mv ' s.Fout 'tmp ' s.Fout]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  rhoMean=buildRhoMeanOperate(threads,s,Dim,f)
    if ~exist(s.Fout,'file')
        initNcFile([s.Fout 'tmp'],'RhoMean',Dim.ws);
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
    end
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
            initNcFile([s.Fout{tt} 'temp'],'density',Dim.ws);
            f.ncvp([s.Fout{tt} 'temp'],'density',f.mDit(RHO{1},Dim.ws),[0 0 0], [Dim.ws]);
            system(['mv ' s.Fout{tt} 'temp '  s.Fout{tt}]);
        end
    end
end
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
