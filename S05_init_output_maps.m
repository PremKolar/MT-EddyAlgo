%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S05_init_output_maps
    %% init
    DD=initialise('');
    %%
    [MAP]=MakeMaps(DD);
    %%
    DD.threads.num=init_threads(DD.threads.num);
    %% find respective index for all grid points of input map    
    MAP.idx=main(DD,MAP);    
    %% save MAP
    save([DD.path.root,'protoMaps.mat'],'-struct','MAP'	)
    %% update infofile
    save_info(DD)
end
function idx=main(DD,MAP)
    if DD.debugmode
        idx=spmd_body(DD,MAP);
    else
        spmd(DD.threads.num)
            idx=spmd_body(DD,MAP) ;
            %% merge composite
            idxx=gop(@vertcat,idx,1);
        end
        idx=sum(idxx{1});
    end
end
function idx=spmd_body(DD,out)
    %% get input example lon/lat
    in.lon=(extractdeepfield(read_fields(DD,1,'cuts'),'grids.LON'));
    in.lat=(extractdeepfield(read_fields(DD,1,'cuts'),'grids.LAT'));
    %% get codisp'ed indeces
    lims=thread_distro(DD.threads.num,numel(in.lon));
    JJ=lims(labindex,1):lims(labindex,2);
    %%
    idx=zeros(1,DD.map.window.size.X*DD.map.window.size.Y);
    %% get Indices For Out Maps
    idx=getIndicesForOutMaps(in,out,JJ,idx);
end
function idx=getIndicesForOutMaps(in,out,JJ,idx)
    %% allocate indices to be calculated by worker
    T=disp_progress('init','allocating old indices to output indeces');
    locSize=numel(JJ);	out.proto=[]; % save mem
    %% loop over indeces
    x=idx;
    y=idx;
    for ii=JJ
        T=disp_progress('disp',T,locSize,100);
        [idx(ii),x(ii),y(ii)]=rangeOp(in.lon(ii),in.lat(ii), out);
    end
    
end
function [lin,x,y]=rangeOp(inLon,inLat,out)
    %% scan for lat/lon within vicinity and use those only
    temp.lon=abs(out.lon-inLon)<=abs(2*(out.inc.x));
    temp.lat=abs(out.lat-inLat)<=abs(2*(out.inc.y));
    used.flag=temp.lon & temp.lat;
    %% out of bounds
    if ~any(used.flag(:)), lin=nan;	return;	end
    %% set lon/lat to be inter-distance checked
    [yi,xi]=find(used.flag);
    used.lat=out.lat(used.flag);
    used.lon=out.lon(used.flag);
    %% find best fit between new/old
    [used.idx]=TransferIdx(inLon,inLat,used);
    %% reset to full size
    y=yi(used.idx.y);
    x=xi(used.idx.x);
    lin=drop_2d_to_1d(y,x,out.dim.y);
end
function [idx]=TransferIdx(lon,lat,used)
    %% build lon/lat matrices
    [yy,xx]=yyxx(lon,lat,used);
    %% take norm2 (sufficient for small distances)
    H=hypot(yy,xx);
    %% find pos of min
    [~,pos]=min(H(:));
    %% raise to 2d to find respective x/y
    [idx.y,idx.x]=raise_1d_to_2d(size(H,1),pos);
end
function [YY,XX]=yyxx(lon,lat,used)
    %% zonal dists
    [A,B]=meshgrid(lon,used.lon);
    xx=abs(A-B)*cosd(lat);
    %% merid dist
    [A,B]=meshgrid(lat,used.lat);
    yy=abs(A-B);
    %%
    [XX,YY]=meshgrid(xx,yy);
end
function [MAP]=MakeMaps(DD)
    %% init output map dim
    xvec=linspace(DD.dim.west,DD.dim.east,DD.dim.X);
    yvec=linspace(DD.dim.south,DD.dim.north,DD.dim.Y);
    [MAP.lon,MAP.lat]=meshgrid(xvec,yvec);
    MAP.proto.nan=nan(size(MAP.lon));
    MAP.proto.zeros=zeros(size(MAP.lon));
    MAP.dim.y=numel(yvec);
    MAP.dim.x=numel(xvec);
    MAP.dim.numel= MAP.dim.y * MAP.dim.x;
    MAP.inc.x=(DD.dim.east-DD.dim.west)/(DD.dim.X-1);
    MAP.inc.y=(DD.dim.north-DD.dim.south)/(DD.dim.Y-1);
end
