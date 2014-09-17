%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S06_init_output_maps
    %% init
    DD=initialise([],mfilename);
    %%
    [MAP]=MakeMaps(DD);
    %%
    DD.threads.num=init_threads(DD.threads.num);
    %% find respective index for all grid points of input map
    MAP.idx=main(DD,MAP);
    %% save MAP
    save([DD.path.root,'protoMaps.mat'],'-struct','MAP'	)
    %% update infofile
    % 	conclude(DD)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx=main(DD,out)
    %% get input example lon/lat
    in.lon=extractdeepfield(read_fields(DD,1,'cuts'),'grids.lon');
    in.lat=extractdeepfield(read_fields(DD,1,'cuts'),'grids.lat');
    idx=zeros(1,DD.map.window.sizePlus.X*DD.map.window.sizePlus.Y);
    spmd
        %% loop over indeces
        T=disp_progress('init','reallocating indeces');
        for ii=drange(1:numel(idx))
            T=disp_progress('disp',T,numel(idx),100);
            [idx(ii)]=getIndicesForOutMaps3(in.lon(ii),in.lat(ii), out);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MAP]=MakeMaps(DD)
    %% init output map dim
    xvec=linspace(DD.map.out.west,DD.map.out.east,DD.map.out.X);
    yvec=linspace(DD.map.out.south,DD.map.out.north,DD.map.out.Y);
    [MAP.lon,MAP.lat]=meshgrid(xvec,yvec);
    MAP.proto.nan=nan(size(MAP.lon));
    MAP.proto.zeros=zeros(size(MAP.lon));
    MAP.dim.y=numel(yvec);
    MAP.dim.x=numel(xvec);
    MAP.dim.numel= MAP.dim.y * MAP.dim.x;
    MAP.inc.x=(DD.map.out.east-DD.map.out.west)/(DD.map.out.X-1);
    MAP.inc.y=(DD.map.out.north-DD.map.out.south)/(DD.map.out.Y-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
