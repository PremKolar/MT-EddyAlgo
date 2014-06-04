%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S06a_getMeanU
    %% init
    DD=initialise;
    %% find files
    [file]=findVelFiles(DD);
    %% get dims
    [d,pos,dim]=getDims(file,DD);
    %% means
    means=getMeans(d,pos,dim,file,DD); %#ok<NASGU>
    %% save
    save([DD.path.meanU.file], 'means')
    disp(['done!'])
end
function means=getMeans(d,pos,dim,file,DD)
    for kk=1:numel(file)
        disp(['found ' file(kk).U ' and ' file(kk).V])
        U(:,:,kk)=squeeze(nc_varget(file(kk).U,'UVEL',dim.start,dim.length))/100; %#ok<*AGROW>
        V(:,:,kk)=squeeze(nc_varget(file(kk).V,'VVEL',dim.start,dim.length))/100;


	x=DD.map.window.size.X;
    y=DD.map.window.size.Y;
    U=downsize(U,x,y);
    V=downsize(V,x,y);

        
    end
    disp(['creating means'])
    U(U<-1e33)=nan; % missing values
    V(V<-1e33)=nan; % missing values
    means.zonal=nanmean(U,3);
    means.merid=nanmean(V,3);
    means.total=hypot(means.zonal,means.merid);
    means.direc=azimuth(zeros(size(means.zonal)),zeros(size(means.zonal)),means.merid,means.zonal);
    means.depth=d(pos.z.start);
    %%
    lin=extractfield(load(DD.path.protoMaps.file),'idx');
    means.small.zonal=nan(DD.map.out.Y,DD.map.out.X);
    for li=unique(lin(lin~=0 & ~isnan(lin)))
        means.small.zonal(li)=nanmean(means.zonal(lin==li));
    end
end
function [d,pos,dim]=getDims(file,DD)
    dWanted=DD.parameters.meanU;
    d=nc_varget(file(1).U,'depth_t');
    [~,pos.z.start]=min(abs(d-dWanted));
    pos.z.start=pos.z.start - 1; % starts at 0
    pos.z.length=1;
    pos.x.start=DD.map.window.limits.west - 1;
    pos.x.length=DD.map.window.size.X;
    pos.y.start=DD.map.window.limits.south-1;
    pos.y.length=DD.map.window.size.Y;
    dim.start = [0 pos.z.start pos.y.start pos.x.start];
    dim.length = 	[inf pos.z.length pos.y.length pos.x.length ];
end
function [file]=findVelFiles(DD)
    %% find the U and V files
    ucc=0; vcc=0;
    for kk=1:numel(DD.path.TempSalt.files)
        if ~isempty(strfind(DD.path.TempSalt.files(kk).name,'UVEL'))
            ucc=ucc+1;
            file(ucc).U=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name]; %#ok<AGROW>
        end
        if ~isempty(strfind(DD.path.TempSalt.files(kk).name,'VVEL'))
            vcc=vcc+1;
            file(vcc).V=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name]; %#ok<AGROW>
        end
    end
end
