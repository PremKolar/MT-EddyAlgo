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
	init_threads(DD.threads.num);
	%% find respective index for all grid points of input map
	spmd(DD.threads.num)
		idx=spmd_body(DD,MAP);
	end
	%% read from composite
	MAP.idx=idx{1};
	%% save MAP
	save([DD.path.root,'protoMaps.mat'],'-struct','MAP'	)
	%% update infofile
	save_info(DD)
end

function idx_out=spmd_body(DD,out)
	%% get input fields
	in.lon=extractdeepfield(read_fields(DD,1,'cuts'),'grids.LON');
	in.lat=extractdeepfield(read_fields(DD,1,'cuts'),'grids.LAT');
	in.idx1d=codistributed(1:numel(in.lon));
	in.dims.y=extractdeepfield(read_fields(DD,1,'cuts'),'window.size.Y');
	%% get Indices For Out Maps
	idx=getIndicesForOutMaps(in,out);
	%% sum vectors from all workers
	idx_out=gop(@horzcat,idx,1);
end


function [MAP]=MakeMaps(DD)
	%% init output map dims
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
	MAP.idx(1,MAP.dim.numel)=struct;
	temp=num2cell(nan(size(MAP.idx)));
	[MAP.idx.x]=deal(temp{1});
	[MAP.idx.y]=deal(temp{1});
	[MAP.idx.lin]=deal(temp{1});
end


function out=getIndicesForOutMaps(in,out)
	%% allocate indices to be calculated by worker
	T=disp_progress('init','allocating old indices to output indeces');
	locSize=numel(getLocalPart(in.idx1d));
	for ii=drange(in.idx1d)
		T=disp_progress('disp',T,locSize,100);
		out=drangeOp(in,out,ii);
	end
end


function out=drangeOp(in,out,ii)
	%% scan for lat/lon within vicinity and use those only
	temp.lon=abs(out.lon-in.lon(ii))<=2*(out.inc.x);
	temp.lat=abs(out.lat-in.lat(ii))<=2*(out.inc.y);
	used.flag=temp.lon & temp.lat;
	%% out of bounds
	if ~any(used.flag),	return;	end
	%% set lon/lat to be inter-distance checked
	[yi,xi]=find(used.flag);
	used.lat=out.lat(used.flag);
	used.lon=out.lon(used.flag);
	%% find best fit between new/old
	[used.idx]=TransferIdx(in.lon(ii),in.lat(ii),used);
	%% reset to full size
	out.idx(ii).y=yi(used.idx.y);
	out.idx(ii).x=xi(used.idx.x);
	out.idx(ii).lin=drop_2d_to_1d(out.idx(ii).y,out.idx(ii).x,in.dims.y);
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
