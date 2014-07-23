%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOW
    %% init
    DD=initialise([],mfilename);
    %% loop over seasons
    seasons={'Summer','Winter'};
    for s = 1:2
        DD.path.full3d.name = DD.path.full3d.(seasons{s}).name;
        minOW.(seasons{s})=main(DD);
    end
    %% save
    save([DD.path.root 'minOW'])
    %% post process
    postProc(DD.path.root,minOW,seasons)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function minOW=main(DD)
    %% set up
    [DD]=maxOWsetUp(DD);
    %% spmd
    metaD=maxOWmain(DD);
    %%
    minOW=maxOWprocess(DD,metaD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function postProc(rootdir,minOW,seasons)
    su=minOW.Summer.ziIntrl;
    wi=minOW.Winter.ziIntrl;
    %% interpolate seasons
    YY=nan([365,size(su)]);
    for ii=1:numel(su)
        YY(1:182,ii)=round(linspace(su(ii),wi(ii),182));
        YY(183:end,ii)=round(linspace(wi(ii),su(ii),183));
    end
    save([rootdir 'ZIfullYear.mat'],'YY')  ;
    %% plot
    for s = 1:2
        plotstuff(minOW.(seasons{s}),minOW.(seasons{1}).depth)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function histDepOw=histUniqDepths(depz,zi,owLg,owAx)    
        histDepOw = histc(owLg( zi == depz ), owAx);
end

function plotstuff(full,depth)
     fuz=reshape(full.zi(~isnan(full.zi)),1,[]);
    uniDepths=unique(fuz);
    %     histc(full.zi(:), uniDepths)  ;
    set(0,'defaulttextinterpreter','latex')
    owLg=-full.ow;
    owAx=[0 logspace(log10(nanmean(owLg(:))),log10(nanmax(owLg(:))),40)];
    histDepOw=nan(numel(uniDepths),numel(owAx));
   
    parfor ud=1:numel(uniDepths)
       histDepOw(ud,:)=histUniqDepths(uniDepths(ud),fuz,owLg,owAx)
    end
    histDepOw(histDepOw==0)=nan;
    bar3(log10(histDepOw));
    view(60.5,28)
    %%
    nyt=ceil(linspace(find(~isnan(depth),1),numel(depth),7));
    nytl =cellfun( @(tk) sprintf('%3.2f',tk), num2cell(depth(nyt)/1000), 'uniformoutput',false)  ;
    set(gca,'ytick',nyt)
    set(gca,'yticklabel',nytl)
    ylabel(['depth [$km$]'])
    %%
    nxt=[find(diff(round(owAx/1e-6)))];
    nxtl=round(owAx(nxt)/1e-6);
    nxtl = cellfun( @(tk) sprintf('-%d',tk), num2cell(nxtl), 'uniformoutput',false)  ;
    xlabel(['$10^6$ Okubo-Weiss Parameter [$1/m^{2}$]'])
    set(gca,'xticklabel',nxtl)
    %%
    zt=get(gca,'ztick')  ;
    nzt=zt;
    nztl=cellfun( @(tk) sprintf('%0.0g',tk), num2cell(nzt), 'uniformoutput',false)  ;
    zlabel(['log10(count)']);
    set(gca,'ztick',nzt)
    set(gca,'zticklabel',nztl);
    saveas(gcf,'yo.fig');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
