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
MAP.idx=getIndicesForOutMaps(DD,CutOne.grids,MAP);
end
function idx_out=getIndicesForOutMaps(DD,grids,MAP)
threads=DD.threads.num;
[Y,X]=size(grids.LAT);
%% copy onto stack for performance
[out.lon]=(MAP.GLO);
[out.lat]=(MAP.GLA);
[out.dlon]=abs(diff(out.lon(1,1:2)));
[out.dlat]=abs(diff(out.lat(1:2)));
%%
in.lon=grids.LON(:);
in.lat=grids.LAT(:);
%%
lims=thread_distro(threads,Y*X);
spmd(DD.threads.num)
%% allocate indices to be calculated by worker
idcs=(lims(labindex,1):lims(labindex,2));
T=disp_progress('init','allocating old indices to output indeces');
kk=0;
for ii=idcs
    kk=kk+1;
    T=disp_progress('disp',T,numel(idcs),10);
    temp.lon=abs(out.lon-in.lon(ii))<=2*(out.dlon);
    temp.lat=abs(out.lat-in.lat(ii))<=2*(out.dlat);
    used.flag=temp.lon & temp.lat;
    if ~any(used.flag) % outofbounds
        idx(kk).y=nan; idx(kk).x=nan; idx(kk).lin=nan;
        continue
    end
    [yi,xi]=find(used.flag);
    used.lat=out.lat(used.flag);
    used.lon=out.lon(used.flag);
    [used.idx]=TransferIdx(in.lon(ii),in.lat(ii),used);
    idx(kk).y=yi(used.idx.y);
    idx(kk).x=xi(used.idx.x);
    idx(kk).lin=drop_2d_to_1d(idx(kk).y,idx(kk).x,MAP.dim.y);
end
%% sum vectors from all workers
idx_out=gop(@horzcat,idx,1);
end
%% send distributed vector to master
idx_out=idx_out{1};

end
function [idx]=TransferIdx(lon,lat,used)

[yy,xx]=yyxx(lon,lat,used);
H=hypot(yy,xx);
[~,pos]=min(H(:));
[idx.y,idx.x]=raise_1d_to_2d(size(H,1),pos) ;

end

function [YY,XX]=yyxx(lon,lat,used)
[A,B]=meshgrid(lon,used.lon);
xx=abs(A-B)*cosd(lat);
[A,B]=meshgrid(lat,used.lat);
yy=abs(A-B);
[XX,YY]=meshgrid(xx,yy);
end
