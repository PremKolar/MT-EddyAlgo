%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD]=maxOWsetUp(DD)
    %% threads
    DD.threads.num=init_threads(DD.threads.num);
    %% find temp and salt files
    [DD.path.TSow]=DataInit(DD);
    %% get window according to user input    
    DD.TSow.window=getfield(load([DD.path.root 'window']),'window');
    %% get z info
    DD.TSow.window=mergeStruct2(DD.TSow.window, GetFields(DD.path.TSow.files(1).salt, cell2struct({DD.TS.keys.depth},'depth')));
    DD.TSow.window.size.Z=numel(DD.TSow.window.depth);
    DD.TSow.window.size.X=DD.map.window.sizePlus.X;
    %%
    DD.path.TSow=appendFields(DD.path.TSow,Data2Init(DD));
    %% distro time steps to threads
    DD.TSow.lims.inTime=thread_distro(DD.threads.num,numel(DD.path.TSow.files));
    DD.TSow.lims.inZ=thread_distro(DD.threads.num,DD.TSow.window.size.Z);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out]=DataInit(DD)
    %% find the temp and salt files
    DirIn=dir([DD.path.full3d.name '*.nc'])   ;
    tt=0;ss=0;
    for kk=1:numel(DirIn);
        if ~isempty(strfind(upper(DirIn(kk).name),DD.TS.keys.salt))
            ss=ss+1;
            out.files(ss).salt=[DD.path.full3d.name DirIn(kk).name];
        end
        if ~isempty(strfind(upper(DirIn(kk).name),DD.TS.keys.temp))
            tt=tt+1;
            out.files(tt).temp=[DD.path.full3d.name DirIn(kk).name];
        end
    end
    out.fnum=numel(out.files);
    out.dir   = [DD.path.Rossby.name];
    out.dailyOWName   = [ out.dir 'OW_'];
    out.dailyRhoName   = [ out.dir 'rho_'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out]=Data2Init(DD)
    TSow=DD.path.TSow;
    dim=DD.TSow.window.size;
    out=initNcRhoMeanPart('densityMean',dim,TSow.dir);
    parfor tt=1:TSow.fnum
        disp(sprintf('init NC file %2d of %2d',tt,TSow.fnum));
        [rho(tt),OW(tt)]=parInitNcs(tt,dim,TSow);
    end
    out.geo=initNcGeoInfo(dim,TSow.dir);
    out.rho=rho;
    out.OW=OW;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho, OW]=parInitNcs(tt,dim,TSow)
    rho={initNcFile('density',dim,tt,TSow.dailyRhoName)};
    OW={initNcFile('OkuboWeiss',dim,tt,TSow.dailyOWName)};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fname=initNcGeoInfo(WinSize,ds)
    fname=[ds  'LatLonDepth.nc' ];
    nc_create_empty(fname,'clobber');
    nc_adddim(fname,'i_index',WinSize.X);
    nc_adddim(fname,'j_index',WinSize.Y);
    nc_adddim(fname,'k_index',WinSize.Z);
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fname=initNcRhoMeanPart(toAdd,WinSize,ds)
    fname.mean=[ds  'RhoMean.nc' ];
    nc_create_empty(fname.mean,'clobber');
    nc_adddim(fname.mean,'k_index',WinSize.Z);
    nc_adddim(fname.mean,'i_index',WinSize.X);
    nc_adddim(fname.mean,'j_index',WinSize.Y);
    %% rho meaned in time
    varstruct.Name = toAdd;
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'k_index','j_index','i_index' };
    nc_addvar(fname.mean,varstruct)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fname=initNcFile(toAdd,WinSize,ff,ds)
    fname=[ds sprintf('%04d.nc',ff) ];
    if ~exist(fname,'file')
        nc_create_empty(fname,'clobber');
        nc_adddim(fname,'k_index',WinSize.Z);
        nc_adddim(fname,'i_index',WinSize.X);
        nc_adddim(fname,'j_index',WinSize.Y);
        %% rho
        varstruct.Name = toAdd;
        varstruct.Nctype = 'double';
        varstruct.Dimension = {'k_index','j_index','i_index' };
        nc_addvar(fname,varstruct)
    end
end
