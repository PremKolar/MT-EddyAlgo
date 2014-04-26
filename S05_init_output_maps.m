%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S05_init_output_maps
	%% init
	DD=initialise('cuts');
	%%
	DD.threads.tracks=thread_distro(DD.threads.num,numel(DD.path.tracks.files));
	DD.CutOne=read_fields(DD,1,'cuts');
	%%
	init_threads(DD.threads.num);
	[MAP]=MakeMaps(DD);
	%% save MAP
	save([DD.path.root,'protoMaps.mat'],'-struct','MAP'	)
	DD.OutMapStuff=MAP;
	%% update infofile
	save_info(DD)
end

function [MAP]=MakeMaps(DD)
	%% init output map dims
	
	xvec=round(10^DD.dim.NumOfDecimals*linspace(DD.dim.west,DD.dim.east,DD.dim.X))/10^DD.dim.NumOfDecimals;
	yvec=round(10^DD.dim.NumOfDecimals*linspace(DD.dim.south,DD.dim.north,DD.dim.Y))/10^DD.dim.NumOfDecimals;
	[MAP.GLO,MAP.GLA]=meshgrid(xvec,yvec);
	MAP.proto.nan=nan(size(MAP.GLO));
	MAP.proto.zeros=zeros(size(MAP.GLO));
	%% find respective index for all grid points of input map
	MAP.idx=getIndicesForOutMaps(DD.CutOne.grids,MAP,DD.threads.num);
end
function idx_out=getIndicesForOutMaps(grids,MAP,threads)
	[Y,X]=size(grids.LAT);
	idx=zeros(Y*X,1);
	%% copy onto stack for performance
	LON=deg2rad(grids.LON);
	LAT=deg2rad(grids.LAT);
	LONout=deg2rad(MAP.GLO);
	LATout=deg2rad(MAP.GLA);
	lims=thread_distro(threads,numel(idx));
	spmd
		%% allocate indices to be calculated by worker
		idcs=(lims(labindex,1):lims(labindex,2));
		T=disp_progress('init','allocating old indices to output indeces');
		for ii=idcs
			T=disp_progress('disp',T,numel(idcs),10);
			idx(ii)=TransferIdx(LON(ii),LAT(ii),LONout,LATout);
		end
		%% sum vectors from all workers
		idx_out=gplus(idx,1);
	end
	
	%% send distributed vector to master
	idx_out=idx_out{1};
end
function [idx]=TransferIdx(lon,lat,LONout,LATout)
	[A,B]=meshgrid(lon,LONout);
	xx=abs(A-B)*cos(lat);
	[A,B]=meshgrid(lat,LATout);
	yy=abs(A-B);
	[~,idx]=min(hypot(xx,yy)) ;
end

