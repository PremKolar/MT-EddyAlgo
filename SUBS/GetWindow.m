function [window,lonlat]=GetWindow(file,mapin,filePattern)
    addpath(genpath('./'));   
    disp('assuming identical LON/LAT for all files!!!')
    %% get data
    % only need lon and lat
    patt.lon=filePattern.lon;
    patt.lat=filePattern.lat;
    [lonlat]=GetFields(file,patt);
    %% find window mask
    window=FindWindowMask(lonlat,mapin);
    %% find rectangle enclosing all applicable data
    window.limits=FindRectangle(window.flag);
    %% size
    window.size=WriteSize(window.limits);
end
function S=WriteSize(lims)
    S.X = lims.east-lims.west   +1;
    S.Y = lims.north-lims.south +1;
end
function window=FindWindowMask(grids,M)
    %% tag all grid points fullfilling all desired lat/lon limits
    if M.east>M.west
        window.flag= grids.lon>=M.west & grids.lon<=M.east & grids.lat>=M.south & grids.lat<=M.north ;
    elseif M.west>M.east  %crossing 180 meridian
        window.flag=((grids.lon>=M.west & grids.lon<=180) | (grids.lon>=-180 & grids.lon<=M.east)) & grids.lat>=M.south & grids.lat<=M.north ;
    end
end
function limits=FindRectangle(flag)
    %% find index limits
    [rows,cols]=find(flag);
    limits.west=min(cols);
    limits.east=max(cols);
    limits.north=max(rows);
    limits.south=min(rows);
end
