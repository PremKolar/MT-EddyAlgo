%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOW
    %% set up
    [DD]=maxOWsetUp;
    %% spmd
    maxOWmain(DD);
    %% make netcdf
    DD=maxOWwriteNCfile(DD);
    %%
    minOW=OWprocess(DD);
    %%
    appendToRossbyNc(DD,minOW);
    %% update DD
    save_info(DD);
end

function appendToRossbyNc(DD,minOW)
    file=DD.path.Rossby.NCfile;
    %% OWmin
    varstruct.Name = 'MinOkuboWeiss';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'j_index','i_index' };
    nc_addvar(file,varstruct)
    %% OWmin depth
    varstruct.Name = 'MinOkuboWeissDepth';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'j_index','i_index' };
    nc_addvar(file,varstruct)
    %%
    [Y,X]=size(minOW.ow);
    nc_varput(file,'MinOkuboWeiss',minOW.ow, [0 0],[Y X]);
    nc_varput(file,'MinOkuboWeissDepth',minOW.z, [0 0],[Y X]);
end