function sub09_scaleStuff
    load S09main II DD T
    velZonmeans(DD,II,T)
    scaleZonmeans(DD,II,T)
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
