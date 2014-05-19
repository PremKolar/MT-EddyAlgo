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
    CutOne=read_fields(DD,1,'cuts');
    %%
    init_threads(DD.threads.num);
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
    [Y,X]=size(grids.LAT);
    idx=zeros(Y*X,1);
    %% copy onto stack for performance
    out.lon=deg2rad(MAP.GLO);
    out.lat=deg2rad(MAP.GLA);
    in.lon=deg2rad(grids.LON);
    in.lat=deg2rad(grids.LAT);
    lims=thread_distro(threads,numel(idx));
%     spmd
        %% allocate indices to be calculated by worker
                idcs=(lims(labindex,1):lims(labindex,2));
        T=disp_progress('init','allocating old indices to output indeces');
        for ii=idcs
            T=disp_progress('disp',T,numel(idcs),10);
            idx(ii)=TransferIdx(in,out,ii);
        end
        %% sum vectors from all workers
        idx_out=gplus(idx,1);
%     end
    %% send distributed vector to master
    idx_out=idx_out{1};
    
end
function [idx]=TransferIdx(in,out,ii)
    
   near=false(size(out.lat));
    FirstNumOfpointsToConsider=min([10,size(out.lat)]);
    temp.grid.lon=abs((in.lon(ii)-out.lon));
    temp.grid.lat=abs((in.lat(ii)-out.lat));
    temp.max.x=max(temp.grid.lon,[],1);
    temp.max.y=max(temp.grid.lat,[],2);    
    demandA=(temp.grid.lon<=temp.max.x(FirstNumOfpointsToConsider));
    demandB=(temp.grid.lat<=temp.max.y(FirstNumOfpointsToConsider));
    near(demandA & demandB)=true;
   [A,B]=meshgridQuick(in.lon(ii),out.lon(near));
    xx=abs(A-B)*cos(in.lat(ii));
    [A,B]=meshgridQuick(in.lat(ii),out.lat(near));
    yy=abs(A-B);  
   H=hypot(yy,xx);
   [~,idxi]=min(H(:)) ;
   near=find(near);
    idx=near(idxi);
end

