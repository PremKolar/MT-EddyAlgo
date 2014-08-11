%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 17-Jul-2014 23:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  OWall=maxOWprocess(DD,metaD)

function yo
    DD=initialise([],mfilename);
    load proc.mat
    OWall=maxOWprocess(DD,metaD)
    gkuzg
end

function  OWall=maxOWprocess(DD,metaD)
    
    try
        system(['rm -r ./critical*']);
    end
    
    NC     =initNC(metaD);
    splits =thread_distro(DD.threads.num,NC.S.Y)-1;
    NC.strt=splits(:,1);
    NC.len =splits(:,2)-NC.strt+1;
    
    
    S=cell2mat(struct2cell(NC.S))';
    NC.OWzi=[DD.path.Rossby.name 'OWzi.nc'];
    NC.OWa=[DD.path.Rossby.name 'OWa.nc'];
    try    %#ok<*TRYNC>
        initNcFile(NC.OWzi,'zi',S([4 2 3]));
    end
    try
        initNcFile(NC.OWa,'OW',S([4 2 3]));
    end
    
    
    spmd
        spmdBcalc(NC);
    end
    %         labBarrier;
    %         switch labindex
    %             case 1
    %                 ALL=spmdBMstrRcv(out,DD.threads.num);
    %             otherwise
    %                 spmdBsendtoMstr(out);
    %         end
    %         labBarrier;
    %     end
    cjghc
    
   
 OWall.zi=nc_varget(NC.OWzi,'zi')  ;
    OWall.ow=nc_varget(NC.OWa,'OW')  ;
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
    
  
    
    
   
  
    
    allOW=OWall.ow;
    depthOW=OWall.depth;
    ziOW=OWall.zi;
    ziIntrp=OWall.ziIntrl;
    ziWeighted=OWall.ziWeighted;
    
    
     save('allOW.mat','allOW','-v7.3')     
     save('zi.mat','ziOW','-v7.3')
 save('ziItnrp.mat','ziIntrp','-v7.3')
     save('ziWeighted.mat','ziWeighted','-v7.3')
    
    
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NC=initNC(metaD)
    NC.geo=metaD.daily.geoOut;
    NC.files(numel(metaD.daily.OWFout)).n=struct;
    for cc=1:numel(metaD.daily.OWFout)
        [dr,fi,ex]=fileparts(metaD.daily.OWFout{cc}) ;
        finishedF{cc}=[dr '/done_' fi ex ];
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
function spmdBcalc(NC)
    TT=NC.S.T;
    in.a=[    0     NC.strt(labindex)   0   ];
    in.b=[   NC.S.Z   NC.len(labindex)    NC.S.X  ];
    mkdirp('./critical')
    mkdirp('./criticalOW')
    for tt=1:TT
        out.a=[    tt-1     NC.strt(labindex)   0   ];
        out.b=[  1   NC.len(labindex)    NC.S.X  ];
        disp(num2str(tt))
        %% get dims for nc read
        
        %% load okubo weiss parameter
        OW=nc_varget(NC.files(tt).n,'OkuboWeiss',in.a,in.b);
        %% kill flags
        OW(abs(OW)>1)=nan;
        %% get bathymetry and dims
        if tt==1
            %         lat=nc_varget(NC.files(tt).n,'lat',in.a(3:4),in.b(3:4));
            %         lon=nc_varget(NC.files(tt).n,'lon',in.a(3:4),in.b(3:4));
            %             lon=nc_varget(NC.files(tt).n,'lon',in.a,in.b);
            %             depth=nc_varget(NC.files(tt).n,'depth');
            %             depth(abs(depth)>1e10)=nan;
            [~,~,~,bath]=getBathym(OW);
            %             out.ow=nan(TT,Yi,X);
            %             out.zi=nan(TT,Yi,X);
        end
        %% min in z
        [owm,zi]=nanmin(OW(2:bath-1,:,:),[], 1);
        zi=zi -1; % correct for (2: ...)
        %% save
        

        
      
        NCcritWrite(NC.OWzi,'zi',zi,out.a,out.b)
        NCcritWrite(NC.OWa,'OW',owm,out.a,out.b)
        
      
        
        
        %         out.ow(tt,:,:)= owm;
        %         out.zi(tt,:,:)= zi;
    end
end

function initNcFile(fname,toAdd,WinSize)
    
    nc_create_empty(fname,'noclobber');
    nc_adddim(fname,'t_index',WinSize(1));
    nc_adddim(fname,'i_index',WinSize(3));
    nc_adddim(fname,'j_index',WinSize(2));
    %% rho
    varstruct.Name = toAdd;
    varstruct.Nctype = 'double';
    varstruct.Dimension = {'t_index','j_index','i_index' };
    nc_addvar(fname,varstruct)
    
    
end









