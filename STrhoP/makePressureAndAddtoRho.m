function makePressureAndAddtoRho
    addpath(genpath('../../netU2/SUBS/'));
    init_threads(12);
    %%
    [lalo,rhoFiles]=getData;
    %%
    [Z,Y,X,DEP,~,G]= inits(lalo);
    %%
    parfor ii=1:numel(rhoFiles)
        opDay(rhoFiles,ii,G,DEP,Y,X,Z)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lalo,rhoFiles]=getData
    rhoFiles=dir('./rho_*.nc');
    lalo.la=nc_varget('LatLonDepth.nc','lat');
    lalo.lo=nc_varget('LatLonDepth.nc','lon');
    lalo.de=nc_varget('LatLonDepth.nc','depth');
end


function  [Z,Y,X,DEP,LAT,G]= inits(lalo)
    [Y,X]=  size(lalo.la);
    [Z]= numel(lalo.de) ;
    DEP=[zeros(1,Y,X); repmat(lalo.de,[1,Y,X])];
    ga=@(M) M(:);
    LAT=(permute(repmat(lalo.la,[1,1,Z]),[3,1,2]));
    G=(reshape(sw_g(LAT(:),ga(DEP(2:end,:,:))),[Z,Y,X]));
    DEP=diff(DEP,1,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opDay(rhoFiles,ii,G,DEP,Y,X,Z)
    Fname=rhoFiles(ii).name;
    data=(nc_varget(Fname,'density'));
    P=cumsum(G.*data.*DEP,1);
    Pname=['pressure_' Fname(5:end)];
    ncOp(Pname,P,X,Y,Z,'pressure',squeeze(DEP(:,1,1)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ncOp(Pname,P,X,Y,Z,fieldname,dep)
    nc_create_empty(Pname,'clobber');
    nc_adddim(Pname,'i_index',X);
    nc_adddim(Pname,'j_index',Y);
    nc_adddim(Pname,'k_index',Z);
    %%
    varstruct.Name = fieldname;
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'k_index','j_index','i_index' };
    nc_addvar(Pname,varstruct);
    %%
    nc_varput(Pname,fieldname,P);
    %%
    varstruct.Name = 'depth';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'k_index'};
    nc_addvar(Pname,varstruct);
    %%
    nc_varput(Pname,'depth',dep);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
