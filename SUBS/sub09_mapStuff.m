function sub09_mapStuff
    load S09main II DD T
    senses=DD.FieldKeys.senses;
    lo=II.lo;
    la=II.la;
    dx=double(extractdeepfield(load([DD.path.cuts.name DD.path.cuts.files(1).name]),'fields.dx'));
    w=load([DD.path.root 'window.mat']);
    X=1:w.window.dimPlus.x;
    Y=1:w.window.dimPlus.y;
    [XX,YY]=meshgrid(X,Y);
    [XXq,YYq]=meshgrid(1:size(lo,2),1:size(lo,1));
    dxq=reshape(griddata(XX(:),YY(:),dx',XXq(:),YYq(:)),size(lo));
    % aviCH=load('../aviN/CC.mat')
    figure
    set(gcf,'visible','off')
    for sense=senses';sen=sense{1};
        %%        
      
        clf
        logFive=@(x) log(x)/log(5);
        VVr=II.maps.(sen).radius.toRo/2;
        VVr(VVr<1e-3)=nan;VVr(VVr>1e3)=nan;
        VV=logFive(VVr);
        pcolor(lo,la,VV);shading flat;colormap([(hsv(8))])
        %         clm=T.radiusToRo;
        clm=[logFive([.125 8]) 9]; % base 5
        decorate(clm,T,DD,sen,'Radius/(2Lr)','km',5,1,1);
        axis(T.axis)   % axis([-180 180 -70 70]);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapRoLLog-' sen],'dpdf');
        %%
        clf
        VV=II.maps.(sen).radius.toRo;
        VV(VV<1e-3)=nan;VV(VV>1e3)=nan;
        pcolor(lo,la,VV);shading flat;colormap(hsv(12))
        clm=[0 6 7];
        decorate(clm,T,DD,sen,'Radius/Lr','km',0,1,1);
        axis(T.axis)   % axis([-180 180 -70 70]);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapRoL-' sen],'dpdf');
        %%
        clf
        VVs=II.maps.(sen).radius.mean.std;
        VVm=II.maps.(sen).radius.mean.mean;
        VV=((VVs./VVm)*100);
        VV(VV<0)=nan;
        VV=log10(VV);
        pcolor(lo,la,VV);shading flat;colormap(hsv(5))
        clm=[log10([1 100 ]) 6];
        decorate(clm,T,DD,sen,' scale: std/mean ','%',10,0,1);
        axis(T.axis)   % axis([-180 180 -70 70]);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapRadStdOMean-' sen],'dpdf');
        %%
        VV=II.maps.(sen).radius.mean.mean/1000;
        colormap(hsv(14));
        pcolor(lo,la,VV);shading flat;colormap(hsv(14));
        clm=[20 160 8];
        decorate(clm,T,DD,sen,'radius','km',0,1,1);
        axis(T.axis)   % axis([-180 180 -70 70]); ;
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapRad-' sen],'dpdf');
        
        %%
        clf
        VV=(II.maps.(sen).radius.mean.mean./dxq);
        a=full(floor(nanmin(VV(:))));
        b=full(ceil(nanmax(VV(:))));
        clm=[a b b-a+1];
        pcolor(lo,la,VV);shading flat;
        colormap(hsv(clm(3)-1));
        %      clm=[20 160 8];
        cb=decorate(clm,T,DD,sen,'radius/dx',' ',0,1,1);
        axis(T.axis)   % axis([-180 180 -70 70]); ;
        xl=(get(cb,'yticklabel'));
        xlc=cell(size(xl,1),1);
        for n=1:size(xl,1);
            xlc{n}=xl(n,:);
            if mod(n,10)~=0
                xlc{n}=' ';
            end
        end
        set(cb,'yticklabel',xlc)
        savefig(DD.path.plots,T.rez,T.width,T.height,['radOdx-' sen],'dpdf');
        %%
        clf
        VV=II.maps.(sen).vel.zonal.mean*100;
        pcolor(lo,la,VV);shading flat
        cw=jet(20);
        cm=[0 0 0];
        ce=(winter(4));
        colormap([cw;cm;ce(:,[1 3 2])])
        decorate([-20 5 6],T,DD,sen,'Zonal velocity','cm/s',0,1,1);
        axis(T.axis)   % axis([-180 180 -70 70]);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapVel-' sen],'dpdf');
        %         CC.(sen).v=VV;
        
        %%
        VV=log(II.maps.(sen).age.mean);
        pcolor(lo,la,VV);shading flat;colormap(jet)
        decorate([log(T.age([1 2])) T.age(3)],T,DD,sen,'age','d',exp(1),0,1);
        axis(T.axis)   % axis([-180 180 -70 70]);
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
        axis(T.axis)   % axis([-180 180 -70 70]);
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
        
        axis(T.axis)   % axis([-180 180 -70 70]);
        savefig(DD.path.plots,T.rez,T.width,T.height,['MapVisitsAll-' sen],'dpdf')
        
        %% TODO
        try
            comp2chelt(II,aviCH)
        end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function comp2chelt(II,aviCH)
    %%
    VV=full(II.maps.(sen).vel.zonal.mean*100);
    VVniko=full(aviCH.CC.(sen).v);
    VVrat=VV-VVniko ;
    LL=VVrat;
    wmax=median(min(LL));
    emax=median(max(LL));
    LL(LL>1 & LL<5)=1.1;
    LL(LL>5 & LL<10)=1.3;
    LL(LL>10 & LL<emax)=1.5;
    LL(LL<-1 & LL>-5)=-1.1;
    LL(LL<-5 & LL>-10)=-1.3;
    LL(LL<-10 & LL>-wmax)=-1.5;
    pcolor(lo,la,LL);shading flat;
    decorate([-log(100) 0 log(100)],T,DD,sen,'U(aviso) [CH-N]','cm/s',0,1,1);
    va=-1.6;vb=1.6;
    cb=colorbar;
    caxis([va vb])
    ct=va:.2:vb;
    ctl=cellfun(@(c) sprintf('%2.1f',c) ,num2cell(ct),'uniformoutput',false);
    ctl(1:3)={'','-10','-5'};
    ctl(end-2:end)={'5','10',''};
    cmaa=flipud([1 0 0; 0 0 1; 0 1 0]);
    cmbb=bone(4);
    CM=[cmaa;(spring(5));flipud(summer(5));flipud(cmbb(1:end-1,:))];
    %     CM=doublemap([va 0 vb],cma(:,[2 1 3]),cmb,[0 0 1],10);
    colormap(CM);
    set(cb,'ytick',ct,'yticklabel',ctl);
    axis(T.axis)   % axis([-180 180 -70 70]); ;
    savefig(DD.path.plots,T.rez,T.width,T.height,['CHmN_aviU-' sen],'dpdf');
    
    %%
    VV=II.maps.(sen).radius.mean.mean/1000;
    VVniko=aviNiko.CC.(sen).L;
    VVdiff=(full(VV-VVniko)./VV)*100 ;
    LL=log(abs(VVdiff)).*sign(VVdiff);
    pcolor(lo,la,LL);shading flat;
    
    decorate([-log(100) 0 log(100)],T,DD,sen,'$\sigma$ [$CH/N$ ratio]','%',0,1,1);
    cb=colorbar;
    caxis([-log(100) log(100)])
    ct=linspace(-log(100),log(100),9);
    ctl=exp(abs(ct)).*sign(ct);
    cma=flipud(jet(50));
    cma(:,1)=.5*cma(:,1) + .5*cma(:,3);
    CM=doublemap([[-log(100) 0 log(100)]],cma(:,[2 1 3]),flipud(jet(50)),[0 0 1],10);
    colormap(CM)
    set(cb,'ytick',ct,'yticklabel',round(ctl))
    axis(T.axis)   % axis([-180 180 -70 70]); ;
    savefig(DD.path.plots,T.rez,T.width,T.height,['CHoN_aviL-' sen],'dpdf');
    
    %%
    %       CC.(sen).L=VV;
    VV=II.maps.(sen).radius.mean.mean/1000;
    VVavi=aviCH.CC.(sen).L;
    VVrat=full(abs(VV./VVavi)) ;
    LL=log(VVrat);
    clf
    pcolor(lo,la,LL);shading flat;
    decorate([-1 0 1],T,DD,sen,'$\sigma$ [pop/aviso ratio]',' ',0,1,1);
    
    axis(T.axis)   % axis([-180 180 -70 70]); ;
    colormap(jet(5));
    cb=colorbar;
    ccc=linspace(log(1/4),log(1),5);
    ccc=[ccc diff(ccc([1 2]))];
    caxis(ccc([1 end]));
    ct=(ccc);
    ctl=rats(exp(ct)',5);
    set(cb,'ytick',ct,'yticklabel',ctl)
    savefig(DD.path.plots,T.rez,T.width,T.height,['POPoAVI_chL-' sen],'dpdf');
    
    
    VV=full(II.maps.(sen).vel.zonal.mean*100);
    VVavi=full(aviCH.CC.(sen).v);
    VVrat=VV-VVavi ;
    LL=VVrat;
    pcolor(lo,la,LL);shading flat;
    decorate([-log(100) 0 log(100)],T,DD,sen,'U(CH) [pop-aviso]','cm/s',0,1,1);
    va=-5;vb=5;
    cb=colorbar;
    caxis([va vb])
    ct=va:1:vb;
    ctl=cellfun(@(c) sprintf('%2.0f',c) ,num2cell(ct),'uniformoutput',false);
    %%
    cma=summer(50);
    cmb=flipud(autumn(50));
    CM=[cma;cmb(:,[ 1 2 2])];
    CM(:,[3])=cos(linspace(-pi*.8,pi*.8,100)').^2;
    colormap(CM)
    set(cb,'ytick',ct,'yticklabel',ctl)
    axis(T.axis)   % axis([-180 180 -70 70]); ;
    savefig(DD.path.plots,T.rez,T.width,T.height,['POPmAVI_chU-' sen],'dpdf');
    
    %%
    %     save CC CC
    
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
