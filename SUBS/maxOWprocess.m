%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 17-Jul-2014 23:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  OWall=maxOWprocess(DD,metaD)

function  maxOWprocess(metaD,DD)
    
    
    
    NC     =initNC(metaD);
    S=cell2mat(struct2cell(NC.S))';
    NC.OWzi=[DD.path.Rossby.name 'OWzi.nc'];
    NC.OWa=[DD.path.Rossby.name 'OWa.nc'];
    initNcFile(NC.OWzi,'zi',S([4 2 3]),{'t_index','j_index','i_index' });
    initNcFile(NC.OWa,'OW',S([4 2 3]),{'t_index','j_index','i_index' });
    
    spmdBcalc(NC,DD.path.Rossby.name);
    
    %
    %     OWall.zi=nc_varget(NC.OWzi,'zi') ;
    %     OWall.ow=nc_varget(NC.OWa,'OW') ;
    %     OWall.depth=nc_varget(NC.geo,'depth');
    %
    %
    %
    %     %% focus on strong neg. okubo weiss
    %     clean=OWall.ow < 5*nanmean(OWall.ow(:));
    %     OWall.ow(~clean) = nan;
    %     OWall.zi(~clean) = nan;
    %     %% ow weighted mean of zi
    %     OWall.owSum      = repmat(nansum(OWall.ow,1),[NC.S.T,1,1]);
    %     OWall.ziWeighted = OWall.ow.*OWall.zi./OWall.owSum;
    %     OWall.meaned.z           = squeeze(nansum(OWall.ziWeighted, 1));
    %     OWall.meaned.ow          = squeeze(nanmean(OWall.ow, 1));
    %     flgd=~squeeze(nansum(clean,1));
    %     [y,x]=size(OWall.meaned.z);
    %     [Xq,Yq]=meshgrid(1:x,1:y);
    %     Xfl=Xq;Yfl=Yq;
    %     Xfl(flgd)=[];
    %     Yfl(flgd)=[];
    %     vq = griddata(Xfl,Yfl,OWall.meaned.z(~flgd),Xq,Yq);
    %     OWall.ziIntrl=round(smooth2a(NeighbourValue(isnan(vq),vq),10));
    %     %    pcolor(vqn);
    %     %    colorbar;
    %
    %     allOW=OWall.ow;
    %     depthOW=OWall.depth;
    %     ziOW=OWall.zi;
    %     ziIntrp=OWall.ziIntrl;
    %     ziWeighted=OWall.ziWeighted;
    %
    %
    %     save('allOW.mat','allOW','-v7.3')
    %     save('zi.mat','ziOW','-v7.3')
    %     save('ziItnrp.mat','ziIntrp','-v7.3')
    %     save('ziWeighted.mat','ziWeighted','-v7.3')
    %
    %
    %
    %
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NC=initNC(metaD)
    NC.geo=metaD.geoOut;
    NC.files(numel(metaD.OWFout)).n=struct;
    for cc=1:numel(metaD.OWFout)
        [dr,fi,ex]=fileparts(metaD.OWFout{cc}) ;
        finishedF{cc}=[dr '/' fi ex ];
    end
    [NC.files(:).n] = deal(finishedF{:});
    
    NC.timeStep = datenum('0815','mmdd');
    NC.S=cell2struct(num2cell(getfield(nc_getvarinfo(NC.files(1).n,'OkuboWeiss'),'Size')),{'Z','Y','X'},2);
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
function daily=initDaily(dirbase,tt,S)
    daily.OWzi=[dirbase sprintf('OWzi_%04d.nc',tt)];
    daily.OWa=[dirbase sprintf('OWa_%04d.nc',tt)];
    initNcFile(daily.OWzi,'zi',S,{'j_index','i_index' });
    initNcFile(daily.OWa,'OW',S,{'j_index','i_index' });
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmdBcalc(NC,dirbase)
    f=funcs;
    %% get bathymetry
    NCf1=NC.files(1).n;
    [~,~,~,bath]=getBathym(nc_varget(NCf1,'OkuboWeiss'));
    %%
    for tt=1:NC.S.T
        daily=initDaily(dirbase,tt,[NC.S.Y NC.S.X]);
        %%
        [zi,owm]=spmdBlck(nc_varget(NC.files(tt).n,'OkuboWeiss'),bath,f);
        nc_varput(daily.OWzi,'zi',squeeze(zi));
        nc_varput(daily.OWa,'OW',squeeze(owm));
        %%
        [zi,owm]=spmdBlck(-nc_varget(NC.files(tt).n,'log10NegOW'),bath,f);
        
        %             OW=-ncvOne(nc_varget(NC.files(tt).n,'log10NegOW'));
        %             %% max in z
        %             [owm,zi]=nanmin(OW(2:bath-1,:,:),[], 1);
        %             zi=gcat(zi-1,3,1); % correct for (2: ...)
        %             owm=gcat(owm,3,1);
        %         end
        %
        
        %         nc_varput(NC.OWzi,'zi',zi,[tt-1 0 0],[1 NC.S.Y NC.S.X])
        %          nc_varput(NC.OWa,'OW',owm,[tt-1 0 0],[1 NC.S.Y NC.S.X])
        %          nc_varput(NC.OWa,'OW',owm)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=funcs
    f.ncvOne = @(A,dim) getLocalPart(codistributed(A,codistributor1d(dim)));
    f.nm = @(OW,bath) nanmin(OW(2:bath-1,:,:),[], 1);
    f.gc = @(a,dim) gcat(squeeze(a),dim,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zi,owm]=spmdBlck(data,bath,f)
    spmd
        %% min in z
        mydata=f.ncvOne(data,3);
        mybath=f.ncvOne(bath,2);
        [owm,zi]=f.nm(mydata,mybath);
        zi=f.gc(zi-1,2); % correct for (2: ...)
        owm=f.gc(owm,2);
    end
    zi=zi{1};
    owm=owm{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(fname,toAdd,WinSize,dimName)
    nc_create_empty(fname,'clobber');
    varstruct.Name = toAdd;
    varstruct.Nctype = 'double';
    for ww=1:numel(WinSize)
        nc_adddim(fname,dimName{ww},WinSize(ww));
    end
    varstruct.Dimension = dimName;
    nc_addvar(fname,varstruct)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








