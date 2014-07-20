%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 17-Jul-2014 23:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  minOW=OWprocess(DD)
    NC=init(DD);
    splits=thread_distro(DD.threads.num,NC.S.X)-1;
    NC.strt=splits(:,1);
    NC.len =splits(:,2)-NC.strt+1;
%     spmd
        out=spmdBcalc(NC);
        labBarrier;
        switch labindex
            case 1
                ALL=spmdBMstrRcv(out,DD.threads.num);
            otherwise
                spmdBsendtoMstr(out);
        end
        labBarrier;
%     end
    minOW=ALL{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NC=init(DD)
    tmp= cellfun(@(f) [DD.path.Rossby.name  f],  extractfield(dir([DD.path.Rossby.name 'OW_*.nc']),'name'),'uniformoutput',false);
    NC.files(numel(tmp)).n=struct;
    [NC.files(:).n] = deal(tmp{:});
    
    NC.S=    cell2struct(num2cell(getfield(nc_getvarinfo(NC.files(1).n,'OkuboWeiss'),'Size')),{'T','Z','Y','X'},2);
    NC.timeStep = nc_varget(NC.files(1).n,'time') +1;
    NC.S.T = numel(NC.files);
    NC.Sname=cell2struct(getfield(nc_getvarinfo(NC.files(1).n,'OkuboWeiss'),'Dimension'),{'T','Z','Y','X'},2);
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
        ALL.z =  [ squeeze(ALL.z),  squeeze(labReceive(tt,10*tt)) ];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmdBsendtoMstr(out)
    labSend(out.ow,1,labindex);
    labSend(out.z, 1,10*labindex);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=spmdBcalc(NC)
    TT=NC.S.T;  
    for tt=1:TT
        %% get dims for nc read
        dim.a=[ 0    0   0     NC.strt(labindex) ];
        dim.b=[ 1    inf inf   NC.len(labindex)  ];
        %% load okubo weiss parameter
        OW=nc_varget(NC.files(tt).n,'OkuboWeiss',dim.a,dim.b);
        %% kill flags
        OW(abs(OW)>1)=nan;
        %% get bathymetry and dims
        if tt==1
        lat=nc_varget(NC.files(tt).n,'lat',dim.a(3:4),dim.b(3:4));
        lon=nc_varget(NC.files(tt).n,'lon',dim.a(3:4),dim.b(3:4));
            %             lon=nc_varget(NC.files(tt).n,'lon',dim.a,dim.b);
            depth=nc_varget(NC.files(tt).n,'depth');
            depth(abs(depth)>1e10)=nan;
            [~,Y,Xi,bath]=getBathym(OW);
            out.ow=nan(TT,Y,Xi);
            out.z=nan(TT,Y,Xi);
        end
        %% min in z
        [owm,zi]=nanmin(OW(2:bath-1,:,:),[], 1);
        zi=zi -1; % correct for (2: ...)
        %% save
        full.ow(tt,:,:)= owm;
        full.zi(tt,:,:) =zi;
    end
    %% focus on strong neg. okubo weiss
    full.ow(full.ow > nanmean(full.ow(:))) = nan;
    %% ow weighted mean of zi
    full.owSum      = repmat(nansum(full.ow,1),[TT,1,1]);
    full.ziWeighted = full.ow.*full.zi./full.owSum;
    out.z           = nansum(full.ziWeighted, 1);
    out.ow          = nanmean(full.ow, 1);
    %% plotstuff
    plotstuff(full,depth)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotstuff(full,depth)
    uniDepths=unique(full.zi(:));
%     histc(full.zi(:), uniDepths)  ;
    set(0,'defaulttextinterpreter','latex')
    owLg=-full.ow;
    owAx=[0 logspace(log10(nanmean(owLg(:))),log10(nanmax(owLg(:))),40)];
    histDepOw=nan(numel(uniDepths),numel(owAx));
    for ud=1:numel(uniDepths)
        depz=uniDepths(ud)   ;
        histDepOw(ud,:)= histc(owLg( full.zi == depz ), owAx);
    end
    histDepOw(histDepOw==0)=nan;
    bar3(log10(histDepOw));
    view(60.5,28)
    xt=get(gca,'xtick');
    yt=get(gca,'ytick');
    zt=get(gca,'ztick')  ;
    nyt=ceil(linspace(find(~isnan(depth),1),numel(depth),7));
    nytl =cellfun( @(tk) sprintf('%3.2f',tk), num2cell(depth(nyt)/1000), 'uniformoutput',false)  ;
    set(gca,'ytick',nyt)
    set(gca,'yticklabel',nytl)
    ylabel(['depth [$km$]'])
    nxt=[find(diff(round(owAx/1e-6)))];
    nxtl=round(owAx(nxt)/1e-6);
    nxtl = cellfun( @(tk) sprintf('-%d',tk), num2cell(nxtl), 'uniformoutput',false)  ;
    xlabel(['$10^6$ Okubo-Weiss Parameter [$1/m^{2}$]'])
    set(gca,'xticklabel',nxtl)
    nzt=zt;
    nztl=cellfun( @(tk) sprintf('%0.0g',tk), num2cell(nzt), 'uniformoutput',false)  ;
    zlabel(['log10(count)']);
    set(gca,'ztick',nzt)
    set(gca,'zticklabel',nztl);
    saveas(gcf,'yo.fig');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















