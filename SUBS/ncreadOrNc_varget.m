function [out] = ncreadOrNc_varget(File,field,dimStart,dimLen)
    if nargin<3
        if exist('nc_varget') %#ok<*EXIST>
            out = nc_varget(File,field);
        elseif exist('ncread')
            out = ncread(File,field);
        else
            error('how to read netcdfs?');
        end
    else
        if exist('nc_varget')
            out = nc_varget(File,field,dimStart,dimLen);
        elseif exist('ncread')
            dimStart = dimStart + 1;
            out = ncread(File,field,dimStart,dimLen);
        else
            error('how to read netcdfs?');
        end
    end
end