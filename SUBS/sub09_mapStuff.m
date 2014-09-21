function sub09_mapStuff
     load S09main II DD T   
 senses=DD.FieldKeys.senses;
 lo=II.lo;
 la=II.la;
 for sense=senses;sen=sense{1};
     
     %% clf
     logFive=@(x) log(x)/log(5);
     VVr=II.maps.(sen).radius.toRo/2;
     VVr(VVr<1e-3)=nan;VVr(VVr>1e3)=nan;
     VV=logFive(VVr);
     pcolor(lo,la,VV);shading flat;colormap([hsv(4); flipud(jet(4))])
     %         clm=T.radiusToRo;
     clm=[logFive([.125 8]) 9]; % base 5
     decorate(clm,T,DD,sen,'Radius/(2Lr)','km',5,1,1);
     axis([-180 180 -70 70])
     savefig(DD.path.plots,T.rez,T.width,T.height,['MapRoLLog-' sen],'dpdf');
     %%
     clf
     VV=II.maps.(sen).radius.toRo;
     VV(VV<1e-3)=nan;VV(VV>1e3)=nan;
     pcolor(lo,la,VV);shading flat;colormap(hsv(12))
     clm=[0 6 7]
     decorate(clm,T,DD,sen,'Radius/Lr','km',0,1,1);
     axis([-180 180 -70 70])
     savefig(DD.path.plots,T.rez,T.width,T.height,['MapRoL-' sen],'dpdf');
     %%
     clf
     VVs=II.maps.(sen).radius.mean.std;
     VVm=II.maps.(sen).radius.mean.mean;
     VV=((VVs./VVm-1)*100);
     VV(VV<0)=nan;
     VV=log10(VV);
     
     pcolor(lo,la,VV);shading flat;colormap(hsv(5))
     clm=[log10([1 100 ]) 6];
     decorate(clm,T,DD,sen,' scale: std/mean ','%',10,0,1);
     axis([-180 180 -70 70])
     savefig(DD.path.plots,T.rez,T.width,T.height,['MapRadStdOMean-' sen],'dpdf');
     %%
     VV=II.maps.(sen).radius.mean.mean/1000;
     pcolor(lo,la,VV);shading flat;colormap(hsv(14))
     
     clm=[20 160 8];
     
     decorate(clm,T,DD,sen,'radius','km',0,1,1);
     axis([-180 180 -70 70])
     
     savefig(DD.path.plots,T.rez,T.width,T.height,['MapRad-' sen],'dpdf');
     %%
     clf
     VV=II.maps.(sen).vel.zonal.mean*100;
     pcolor(lo,la,VV);shading flat
     doublemap([T.vel(1),0,T.vel(2)]	,II.aut,II.win,[.9 1 .9],20);
     decorate(T.vel,T,DD,sen,'Zonal velocity','cm/s',0,1,1);
     axis([-180 180 -70 70])
     savefig(DD.path.plots,T.rez,T.width,T.height,['MapVel-' sen],'dpdf');
     %%
     VV=log(II.maps.(sen).age.mean);
     pcolor(lo,la,VV);shading flat;colormap(jet)
     decorate([log(T.age([1 2])) T.age(3)],T,DD,sen,'age','d',exp(1),0,1);
     axis([-180 180 -70 70])
     savefig(DD.path.plots,T.rez,T.width,T.height,['MapAge-' sen],'dpdf')
     %%
     VV=(II.maps.(sen).visits.single);
     VV(VV==0)=nan;
     pcolor(lo,la,VV);shading flat;colormap(jet(11))
     cb=decorate([0 27.5,11],T,DD,sen,'Visits of unique eddy',' ',0,1,1);
%      decorate(T.visitsunique,T,DD,sen,'Visits of unique eddy',' ',0,1,1);
    set(cb,'ytick',[0 5:5:25]) 
    set(cb,'yticklabel',[1 5:5:25]) 
    set(cb,'ylim',[0 27.5]) 
axis([-180 180 -70 70])
savefig(DD.path.plots,T.rez,T.width,T.height,['MapVisitsUnique-' sen],'dpdf')
     %%
    clf
    VV=(II.maps.(sen).visits.all);
     VV(VV==0)=nan;
     pcolor(lo,la,VV);shading flat;colormap(hsv(21))
    cb=decorate([0 105,11],T,DD,sen,'Total Visits',' ',0,1,1);
  set(cb,'ytick',[0 10:10:100]) 
    set(cb,'yticklabel',[1 10:10:100]) 
    set(cb,'ylim',[0 105]) 
     
     axis([-180 180 -70 70])
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawcoast
    load coast;
    hold on; plot(long,lat,'LineWidth',0.5);
end
