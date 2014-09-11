%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Sep-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S09_plotsNew
   DD=initialise([],mfilename);
   save DD
%     clf
%     clear all
    %     initialise
%     dbstop if error
%     load DD
    ticks.rez=get(0,'ScreenPixelsPerInch');
    ticks.width=800;
    ticks.height=600;
    ticks.y= 0;
    ticks.x= 0;
    ticks.age=[1,2*365,10];
    %     ticks.isoper=[DD.thresh.shape.iq,1,10];
    ticks.isoper=[.6,1,10];
    ticks.radius=[50,250,11];
    ticks.radiusStd=[0,150,11];
    ticks.radiusToRo=[1,5,5];
    ticks.amp=[1,20,7];
    %ticks.visits=[0,max([maps.AntiCycs.visitsSingleEddy(:); maps.Cycs.visitsSingleEddy(:)]),5];
    ticks.visits=[1,20,11];
    ticks.visitsunique=[1,10,10];
    ticks.dist=[-1500;500;11];
    %ticks.dist=[-100;50;16];
    ticks.disttot=[1;3000;5];
    ticks.vel=[-30;20;6];
    ticks.axis=[DD.map.out.west DD.map.out.east DD.map.out.south DD.map.out.north];
    ticks.lat=[ticks.axis(3:4),5];
    ticks.minMax=cell2mat(extractfield( load([DD.path.analyzed.name, 'vecs.mat']), 'minMax'));
    %%
    main(DD,ticks)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,T)
    II=initStuff(DD);
    %% collect tracks
    trackstuff(DD,T)
    %%
    %     velZonmeans(DD,II,T)
    %     scaleZonmeans(DD,II,T)
    %     %%
    %     mapStuff(DD,T,II)
end


function trackstuff(DD,T)    
    tsenses=fieldnames(DD.path.analyzedTracks)';
    senses=DD.FieldKeys.senses;
    for ss=1:2
       clear single cats
        tsen=tsenses{ss};
        sen = senses{ss};
        root=DD.path.analyzedTracks.(tsen).name;
        eds= DD.path.analyzedTracks.(tsen).files;
      parfor ff=1:numel(eds)          
          single(ff)=load([root eds(ff).name]);
        end
        for fn=fieldnames(single)'; fn=fn{1};
            cats.(fn)=extractdeepfield(single,fn);
        end
        TR.(sen).single=single;
        TR.(sen).cats=cats;
    end
    %%
    
            jjjjjd
    
     for ss=1:2
    
         tsen=tsenses{ss};
        sen = senses{ss};
       Rm=TR.(sen).cats.radiusmean
       Rm=TR.(sen).cats.median
    end
    
    
    
end

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function velZonmeans(DD,II,T)
    %%
    pp(1)=plot(II.la(:,1),II.maps.zonMean.Rossby.small.phaseSpeed	); 	hold on
    acv=squeeze(nanmean(II.maps.AntiCycs.vel.zonal.mean,2));
    cv=squeeze(nanmean(II.maps.Cycs.vel.zonal.mean,2));
    pp(2)=plot(II.la(:,1),acv	,'r');
    pp(3)=plot(II.la(:,1),cv,'black');
    axis([DD.map.out.south DD.map.out.north min([min(acv) min(cv)]) max([max(acv) max(cv)]) ])
    set(pp(:),'linewidth',2)
    leg=legend('Rossby-wave phase-speed',2,'anti-cyclones net zonal velocity',2,'cyclones net zonal velocity');
    set( get(leg,'children'),'linewidth',2)
    ylabel('[cm/s]')
    xlabel('[latitude]')
    title(['velocity - zonal means [cm/s]'])
    set(get(gcf,'children'),'linewidth',2)
    savefig(DD.path.plots,T.rez,T.width,T.height,['S-velZonmean'],'dpdf');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scaleZonmeans(DD,II,T)
    %%
    pp(1)=plot(II.la(:,1),2*II.maps.zonMean.Rossby.small.radius); 	hold on
    pp(2)=plot(II.la(:,1),II.maps.zonMean.AntiCycs.radius.mean.mean,'r');
    pp(3)=plot(II.la(:,1),II.maps.zonMean.Cycs.radius.mean.mean,'black');
    set(pp(:),'linewidth',2)
    leg=legend('2 x Rossby Radius','anti-cyclones radius','cyclones radius');
    set( get(leg,'children'),'linewidth',2)
    ylabel('[m]')
    xlabel('[latitude]')
    title(['scale - zonal means'])
    maxr=nanmax(nanmax([II.maps.zonMean.AntiCycs.radius.mean.mean(:) II.maps.zonMean.Cycs.radius.mean.mean(:)]));
    axis([min(II.la(:,1)) max(II.la(:,1)) 0 maxr])
    savefig(DD.path.plots,T.rez,T.width,T.height,['S-scaleZonmean'],'dpdf');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function II=initStuff(DD)
    II.aut=autumn(100);
    II.win=winter(100);
    II.maps=load([DD.path.analyzed.name, 'maps.mat']);
    II.la=II.maps.Cycs.lat;
    II.lo=II.maps.Cycs.lon;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mapStuff(DD,T,II)
    senses=DD.FieldKeys.senses;
    for sense=senses;sen=sense{1};
        %%
        VVr=II.maps.(sen).radius.toRo;
        VVr(VVr<1e-3)=nan;VVr(VVr>1e3)=nan;
        VV=log(VVr)/log(5);
        pcolor(lo,la,VV);shading flat;colormap(jet)
        clm=T.radiusToRo;
        clm([ 1 2 ])=log(clm([ 1 2 ]))./log([5 5]); % base 5
        decorate(clm,T,DD,sen,'Radius/Lr','km',5,1,1);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapRoL-' sen],'dpdf');
        %%
        VV=II.maps.(sen).radius.mean.std/1000;
        pcolor(lo,la,VV);shading flat;colormap(jet)
        decorate('radiusStd',T,DD,sen,'Radius std','km',0,1);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapRadStd-' sen],'dpdf');
        %%
        VV=II.maps.(sen).radius.mean.mean/1000;
        pcolor(lo,la,VV);shading flat;colormap(jet)
        decorate('radius',T,DD,sen,'Radius','km',0,1);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapRad-' sen],'dpdf');
        %%
        VV=II.maps.(sen).vel.zonal.mean*100;
        pcolor(lo,la,VV);shading flat
        doublemap([T.vel(1),0,T.vel(2)]	,II.aut,II.win,[.9 1 .9],20);
        decorate('vel',T,DD,sen,'Zonal velocity','cm/s',0,1);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapVel-' sen],'dpdf');
        %%
        VV=log(II.maps.(sen).age.mean);
        pcolor(lo,la,VV);shading flat;colormap(jet)
        decorate('age',T,DD,sen,'age','d',1,1);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapAge-' sen],'dpdf')
        %%
        VV=(II.maps.(sen).visits.single);
        VV(VV==0)=nan;
        pcolor(lo,la,VV);shading flat;colormap(jet)
        decorate('visitsunique',T,DD,sen,'Visits of unique eddy',' ',0,1);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapVisitsUnique-' sen],'dpdf')
        %%
        VV=(II.maps.(sen).visits.all);
        VV(VV==0)=nan;
        pcolor(lo,la,VV);shading flat;colormap(jet)
        decorate('visits',T,DD,sen,'Total Visits',' ',0,1);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapVisitsAll-' sen],'dpdf')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb=decorate(clm,ticks,DD,tit,tit2,unit,logbase,decim,coast)
    %%
    dec=10.^decim;
    %%
    axis(ticks.axis);
    set(gca,'ytick',ticks.y);
    set(gca,'xtick',ticks.x);
    cb=colorbar;
    %%
    zticks=linspace(clm(1),clm(2),clm(3))';
    %%
    switch logbase
        case 0
            zticklabel=num2str(round(zticks*dec)/dec);
        otherwise
            ztl=logbase.^zticks;
            [zaehler,nenner]=rat(ztl);
            nenn=nenner(1);
            s=zaehler>=nenner;
            ztlA=round(10*zaehler(~s).*repmat(nenn,size(nenner(~s)))./nenner(~s))/10;
            zticklabelA=cellfun(@(c) [num2str(c),'/',num2str(nenn)], num2cell(ztlA),'uniformoutput',false);
            ztlB=round(dec.*zaehler(s)./nenner(s))/dec;
            zticklabelB=cellfun(@(c) num2str(c),num2cell(ztlB),'uniformoutput',false);
            zticklabel=[zticklabelA;zticklabelB];
    end
    %%
    caxis([zticks(1) zticks(end)])
    set(cb,'ytick',zticks);
    set(cb,'yticklabel',zticklabel);
    title([tit,' - ',tit2,' [',unit,']'])
    xlabel(['Eddies that died younger ',num2str(DD.thresh.life),' days are excluded'])
    if coast
        drawcoast;
    end
end
%------------------
function [zticklabel,zticks]=logcase(ticks,ratsticks,decim)
    zticks=logspace(log10(ticks(1)),log10(ticks(2)),ticks(3))';
    if ratsticks
        zticklabel=ratscase(exp(zticks));
    else
        zticklabel=round(exp(zticks)*decim)/decim;
        zticklabel=num2str(zticklabel);
    end
end
%----------------
function zticklabel=ratscase(zticks)
    [n,d]=rat(zticks);
    nc=(num2cell(n));
    dc=(num2cell(d));
    zticklabel=cellfun(@(a,b) [num2str(a) '/' num2str(b)],nc,dc,'uniformoutput',false);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawcoast
    load coast;
    hold on; plot(long,lat,'LineWidth',0.5);
end


























