function sub09_scaleStuff
    load S09main II DD T
    velZonmeans(DD,II,T)
    scaleZonmeans(DD,II,T)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function velZonmeans(DD,II,T)
    %%
    clf
    lw=2;
    pp(1)=plot(II.la(:,1),-II.maps.zonMean.Rossby.small.phaseSpeed*100); 	hold on
  
   acv=-nansum(II.maps.AntiCycs.vel.zonal.mean.*II.maps.AntiCycs.visits.all,2)./nansum(II.maps.AntiCycs.visits.all,2)*100
   cv=-nansum(II.maps.Cycs.vel.zonal.mean.*II.maps.Cycs.visits.all,2)./nansum(II.maps.Cycs.visits.all,2)*100
    
    acv=-squeeze(nanmean(II.maps.AntiCycs.vel.zonal.mean,2))*100;
    cv=-squeeze(nanmean(II.maps.Cycs.vel.zonal.mean,2))*100;
    pp(2)=plot(II.la(:,1),acv	,'r');
    pp(3)=plot(II.la(:,1),cv,'black');
    pp(4)=plot([DD.map.out.south DD.map.out.north], [0 0],'b--')
%     axis([DD.map.out.south DD.map.out.north min([min(acv) min(cv)]) max([max(acv) max(cv)]) ])
    axis([DD.map.out.south DD.map.out.north -5 10])
    set(pp(:),'linewidth',lw)
    leg=legend('Rossby-wave phase-speed',2,'anti-cyclones net zonal velocity',2,'cyclones net zonal velocity');
    set( get(leg,'children'),'linewidth',lw)
    ylabel('[cm/s]')
    xlabel('[latitude]')
    title(['westward propagation [cm/s]'])
    set(get(gcf,'children'),'linewidth',lw)
    savefig(DD.path.plots,T.rez,T.width,T.height,['S-velZonmean'],'dpdf');
    
    
    clf  
    hold on  
    lw=2;
    acm=II.maps.AntiCycs.vel.zonal.mean;
    aci=II.maps.AntiCycs.visits.all;
   ccm=II.maps.Cycs.vel.zonal.mean;
    cci=II.maps.Cycs.visits.all;
  CCAV=-nanmean( ( acm .* aci + ccm .* cci ) ./ (aci+cci) ,2)*100;
   
    pp=plot(II.la(:,1),CCAV,'--','color',[0 .5 0]);
    axis([-50 50 -5 20]);
     set(gca,'YAxisLocation','right','ytick',[-5 0 5 10 15 20],'xtick',-50:25:50);
   set(gca,'xAxisLocation','top')
     set(pp,'linewidth',lw)   
%     ylabel('[cm/s]')
%     xlabel('[latitude]')
%     title(['westward propagation [cm/s]'])
    set(get(gcf,'children'),'linewidth',lw)
    grid on
 
    
    savefig(DD.path.plots,100,500,500,['S-velZonmean4chelt11comp'],'dpdf'); 
    
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scaleZonmeans(DD,II,T)
    %% 
    lw=1.5;
    pp(1)=plot(II.la(:,1),2*II.maps.zonMean.Rossby.small.radius); 	hold on
    pp(2)=plot(II.la(:,1),II.maps.zonMean.AntiCycs.radius.mean.mean,'r');
    pp(3)=plot(II.la(:,1),II.maps.zonMean.Cycs.radius.mean.mean,'black');
    set(pp(:),'linewidth',lw)
    leg=legend('2 x Rossby Radius','anti-cyclones radius','cyclones radius');
    set( get(leg,'children'),'linewidth',lw)
    ylabel('[m]')
    xlabel('[latitude]')
    title(['scale - zonal means'])
    maxr=nanmax(nanmax([II.maps.zonMean.AntiCycs.radius.mean.mean(:) II.maps.zonMean.Cycs.radius.mean.mean(:)]));
    axis([min(II.la(:,1)) max(II.la(:,1)) 0 maxr])
    savefig(DD.path.plots,T.rez,T.width,T.height,['S-scaleZonmean'],'dpdf');
end
