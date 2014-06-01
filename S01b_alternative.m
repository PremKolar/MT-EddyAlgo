%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 31-Mai-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get rossby stuff from POP and give to AVISO
function S01b_alternative
    %% init
    DD=initialise('conts');
    %% read input from POP data
    [RS]=getRossbyStuff(DD);
    %% create reshaped data
    saveData(RS)
end


function [RS]=getRossbyStuff(DD)
    RS.file=[DD.path.Rossby.name DD.path.Rossby.files.name];
    RS.c=nc_varget(RS.file,'RossbyPhaseSpeed');
    RS.Lr=nc_varget(RS.file,'RossbyRadius');
    x=DD.map.window.size.X;
    y=DD.map.window.size.Y;
    RS.c=downsize(RS.c,x,y);
    RS.Lr=downsize(RS.Lr,x,y);
end


function saveData(RS)
    %% init
    system(['cp ' RS.file ' ' RS.file '_temp'])
    nc_create_empty(RS.file,'clobber')
    [Y,X]=size(RS.c);
    dim.start  = [0 0];
    dim.length = [Y X];
    nc_adddim(RS.file,'i_index',X);
    nc_adddim(RS.file,'j_index',Y);
    %% Ro1
    varstruct.Name = 'RossbyRadius';
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'j_index','i_index' };
    nc_addvar(RS.file,varstruct)
    nc_varput(RS.file,'RossbyRadius',RS.Lr,dim.start, dim.length);
    %% c1
    varstruct.Name = 'RossbyPhaseSpeed';
    varstruct.Nctype = 'double';
    varstruct.Dimension ={'j_index','i_index' };
    nc_addvar(RS.file,varstruct)
    nc_varput(RS.file,'RossbyPhaseSpeed',RS.c,dim.start, dim.length);
end