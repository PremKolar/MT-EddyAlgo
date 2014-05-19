%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S05_init_output_maps
	%% init
	DD=initialise('');
	DD=load([DD.path.root,'DD.mat'])
	%%
	DD.threads.tracks=thread_distro(DD.threads.num,numel(DD.path.tracks.files));
	CutOne=read_fields(DD,1,'cuts');
	%%
	%     init_threads(DD.threads.num);
	[MAP]=MakeMaps(DD,CutOne); %#ok<*NASGU>
	%% save MAP
	save([DD.path.root,'protoMaps.mat'],'-struct','MAP'	)
	%% update infofile
	save_info(DD)
end

function [MAP]=MakeMaps(DD,CutOne)
	%% init output map dims
	xvec=round(10^DD.dim.NumOfDecimals*linspace(DD.dim.west,DD.dim.east,DD.dim.X))/10^DD.dim.NumOfDecimals;
	yvec=round(10^DD.dim.NumOfDecimals*linspace(DD.dim.south,DD.dim.north,DD.dim.Y))/10^DD.dim.NumOfDecimals;
	[MAP.GLO,MAP.GLA]=meshgrid(xvec,yvec);
	MAP.proto.nan=nan(size(MAP.GLO));
	MAP.proto.zeros=zeros(size(MAP.GLO));
	MAP.dim.y=numel(yvec);
	MAP.dim.x=numel(xvec);
	MAP.dim.numel= MAP.dim.y * MAP.dim.x;
	%% find respective index for all grid points of input map
	MAP.idx=getIndicesForOutMaps(CutOne.grids,MAP,DD.threads.num);
end
function idx_out=getIndicesForOutMaps(grids,MAP,threads)
	[Yout,Xout]=size(MAP.GLO);
	[Y,X]=size(grids.LAT);
	%% copy onto stack for performance
	[out.lon]=deg2rad(MAP.GLO);
	[out.lat]=deg2rad(MAP.GLA);
	%%
	temp=num2cell(deg2rad(grids.LON));
	[in(1:Y*X).lon]=deal(temp{:});
	temp=num2cell(deg2rad(grids.LAT));
	[in(1:Y*X).lat]=deal(temp{:});
	%%
	lims=thread_distro(threads,Y*X);
	%          spmd
	%% allocate indices to be calculated by worker
	idcs=(lims(labindex,1):lims(labindex,2));
	idx(numel(idcs)).lin=nan;idx(numel(idcs)).x=nan;idx(numel(idcs)).y=nan;
	T=disp_progress('init','allocating old indices to output indeces');
	for ii=idcs
		kk=ii-idcs(1)+1;
		T=disp_progress('disp',T,numel(idcs),10);
		[idx(kk).lin,idx(kk).y,idx(kk).x]=TransferIdx(in(ii),out,Yout,Xout);
	end
	%% sum vectors from all workers
	idx_out=gop(@horzcat,idx,1);
	%             end
	%% send distributed vector to master
	idx_out=idx_out{1};
	
end
function [lin,y,x]=TransferIdx(in,OUT,Y,X)
	
	[yy,xx]=yyxx(in,OUT)
	
	
	H=hypot(YY(:),XX(:));
	[Hs]=min(reshape(H,[Y X Y*X]),[],3);
	[Hs]= reshape(H,[Y X Y*X])
	[Hs]=min(reshape(H,[Y Y*X X]),[],2);
	[~,lin]=min(Hs(:))   ;
	[y,x]=raise_1d_to_2d(Y,lin) ;
end

function [yy,xx]=yyxx(in,OUT)
	[A,B]=meshgrid(in.lon,OUT.lon);
	xx=abs(A-B)*cos(in.lat);
	[A,B]=meshgrid(in.lat,OUT.lat);
	yy=abs(A-B);
	%  [YY,XX]=meshgrid(yy,xx);
	
end
