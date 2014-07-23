%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 17-Jul-2014 23:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  OWall=maxOWprocess(DD,metaD)
    NC=init(metaD);
    splits=thread_distro(DD.threads.num,NC.S.Y)-1;
    NC.strt=splits(:,1);
    NC.len =splits(:,2)-NC.strt+1;
    spmd
        out=spmdBcalc(NC);
        labBarrier;
        switch labindex
            case 1
                ALL=spmdBMstrRcv(out,DD.threads.num);
            otherwise
                spmdBsendtoMstr(out);
        end
        labBarrier;
    end
    OWall=ALL{1};
    OWall.depth=nc_varget(NC.geo,'depth');
    %% focus on strong neg. okubo weiss
    clean=OWall.ow < 5*nanmean(OWall.ow(:));
    OWall.ow(~clean) = nan;
    OWall.zi(~clean) = nan;
    %% ow weighted mean of zi
    OWall.owSum      = repmat(nansum(OWall.ow,1),[NC.S.T,1,1]);
    OWall.ziWeighted = OWall.ow.*OWall.zi./OWall.owSum;
    OWall.meaned.z           = squeeze(nansum(OWall.ziWeighted, 1));
    OWall.meaned.ow          = squeeze(nanmean(OWall.ow, 1));
    flgd=~squeeze(nansum(clean,1));
    [y,x]=size(OWall.meaned.z);
    [Xq,Yq]=meshgrid(1:x,1:y);
    Xfl=Xq;Yfl=Yq;
    Xfl(flgd)=[];
    Yfl(flgd)=[];
    vq = griddata(Xfl,Yfl,OWall.meaned.z(~flgd),Xq,Yq);
    OWall.ziIntrl=round(smooth2a(NeighbourValue(isnan(vq),vq),10));
    %    pcolor(vqn);
    %    colorbar;
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NC=init(metaD)
    NC.geo=metaD.daily.geoOut;
    NC.files(numel(metaD.daily.OWFout)).n=struct;
    [NC.files(:).n] = deal(metaD.daily.OWFout{:});
    NC.S=  metaD.dim.ws;
    NC.timeStep = datenum('0815','mmdd');
    NC.S.T = numel(NC.files);
    NC.Sname=cell2struct(getfield(nc_getvarinfo(NC.files(1).n,'OkuboWeiss'),'Dimension'),{'Z','Y','X'},2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z,Y,X,bath]=getBathym(OW)
    [Z,Y,X]=size(OW);
    OW2d=reshape(OW,[Z,Y*X]);
    [~,bathUpdown]=min(isnan(flipud(OW2d)),[],1);
    bath=reshape( Z-bathUpdown + 1, [Y,X]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ALL=spmdBMstrRcv(ALL,threads)
    for tt=2:threads
        ALL.ow = [ squeeze(ALL.ow), squeeze(labReceive(tt,tt)) ];
        ALL.zi =  [ squeeze(ALL.zi),  squeeze(labReceive(tt,10*tt)) ];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmdBsendtoMstr(out)
    labSend(out.ow,1,labindex);
    labSend(out.zi, 1,10*labindex);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=spmdBcalc(NC)
    TT=NC.S.T;
    for tt=1:TT
        %% get dims for nc read
        dim.a=[    0     NC.strt(labindex)   0   ];
        dim.b=[   inf   NC.len(labindex)    inf  ];
        %% load okubo weiss parameter
        OW=nc_varget(NC.files(tt).n,'OkuboWeiss',dim.a,dim.b);
        %% kill flags
        OW(abs(OW)>1)=nan;
        %% get bathymetry and dims
        if tt==1
            %         lat=nc_varget(NC.files(tt).n,'lat',dim.a(3:4),dim.b(3:4));
            %         lon=nc_varget(NC.files(tt).n,'lon',dim.a(3:4),dim.b(3:4));
            %             lon=nc_varget(NC.files(tt).n,'lon',dim.a,dim.b);
            %             depth=nc_varget(NC.files(tt).n,'depth');
            %             depth(abs(depth)>1e10)=nan;
            [~,Yi,X,bath]=getBathym(OW);
            out.ow=nan(TT,Yi,X);
            out.zi=nan(TT,Yi,X);
        end
        %% min in z
        [owm,zi]=nanmin(OW(2:bath-1,:,:),[], 1);
        zi=zi -1; % correct for (2: ...)
        %% save
        out.ow(tt,:,:)= owm;
        out.zi(tt,:,:)= zi;
    end
end













