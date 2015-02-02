function sub09_TPzStuff
    load S09main DD T
    %     [procData]=inits(DD);
    %     save([DD.path.analyzed.name 'procData.mat'],'procData');
    load([DD.path.analyzed.name 'procData.mat'],'procData');
    senses = DD.FieldKeys.senses';
    %     spmd(2)
    %         if labindex==1
    %             %             TPz(DD,T,procData.tracks,senses,'lat',100,'lat',0);
%     TPz(DD,T,procData.tracks,senses,'age',500,'age',1);
    %         else
    TPzGlobe(DD,T,procData.tracks,senses,'age',365,'age',1,1);
    %         end
    %     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TPzGlobe(DD,ticks,tracks,senses,colorfield,minlen,cticks,logornot,fac)
%     close all
    globe=true;
    %     drawLinez(tracks.(sen),minlen)
    cmap = winter(100);
    drawColorLinez(ticks,tracks.(senses{1}),colorfield,minlen,cticks,logornot,globe,fac,cmap) ;
    cb{1} = colorbar;
    colormap(cb{1},winter(100));
    hold on
    %%
    cmap=autumn(100);
    drawColorLinez(ticks,tracks.(senses{2}),colorfield,minlen,cticks,logornot,globe,fac,cmap) ;
    cb{2}=colorbar('westoutside');
    colormap(cb{2},autumn(100))
    
    %%
    axis([-180 180 -70 70])
    drawcoast
    
    tit=['tracks-' colorfield];
    if logornot
        ticks.(cticks)(1:2) = log(ticks.(cticks)(1:2));
        ticks.(cticks)(ticks.(cticks)==0) = 1;
    end
%     set(cb{1},'ytick',0:.2:.8)
    for ss=1:2
        set(cb{ss},'ytick',linspace(0,1,ticks.(cticks)(3)))
        if logornot
            set(cb{ss},'yticklabel',round(exp(linspace(ticks.(cticks)(1),ticks.(cticks)(2),ticks.(cticks)(3)))))
        end
    end
    set(cb{2},'yticklabel',[])
    
    grid on
    
%     titalt=[tit '-alt'];
%     savefig(DD.path.plots,ticks.rez,1400,600,titalt,'dpdf')
    
    dims=[8*12 6*12];FS=20;
    set(gcf,'Visible','off','position',[0 0 1200 800],'paperunits','inch','papersize',dims,'paperposition',[0 0 dims]);
    set(findall(gcf,'type','text'),'FontSize',FS,'interpreter','latex')
    set(gca,'FontSize',FS);
    fname = [DD.path.plots tit];
    print('-depsc','-painters',fname)
    print('-depsc','-painters',fname)
    system(['epstopdf ' fname '.eps']);
    system(['pdfcrop --margins "1 1 1 1" ' fname '.pdf ' fname '.pdf']);
    cpPdfTotexMT(tit);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TPz(DD,ticks,tracks,senses,colorfield,minlen,cticks,logornot)
    close all
    set(gca,'yaxisLocation','right');
    axis([-5000 3000 -2000 2000]);
    %     axis equal
    set(gca,'ytick',[-1000   0   1000]);
    set(gca,'xtick',[-4000 -2000 -1000  0 1000]);
    set(gca,'yticklabel',[-1   0   1]);
    set(gca,'xticklabel',[-4 -2 -1 0 1 ]);
    tit=['defl-' colorfield];
    title(['cyclones in red / a.-cyc. in blue [days]']);
    xlabel('[Mm]');
    grid on;
    cbticks=10;
    %%
    cmap = flipud(parula(100));
    [minV, maxV]=drawColorLinez(ticks,tracks.(senses{1}),colorfield,minlen,cticks,logornot,0,1,cmap) ;
    cb{1} = colorbar;
    set(cb{1},'ytick',linspace(0,1,cbticks))
    set(cb{1},'yticklabel',linspace(minV,maxV,cbticks))
    colormap(cb{1},cmap);
    hold on
    
    cmap=flipud(hot(110));
    cmap=cmap(11:end,:);
    [minV, maxV]=drawColorLinez(ticks,tracks.(senses{2}),colorfield,minlen,cticks,logornot,0,1,cmap) ;
    cb{2}=colorbar('westoutside');
    p2 = get(cb{2},'position')
    p2(3) = p2(3)/2;
    p2(1) = p2(1) -  p2(3);
    set(cb{2},'position',p2);
    p1 = get(cb{1},'position');
    p1(3) = p1(3)/2;
    p1(1) = p2(1) -  p2(3);
    set(cb{1},'position',p1);
    set(cb{2},'ytick',linspace(0,1,cbticks));
    set(cb{2},'yticklabel','');
    colormap(cb{2},cmap);
    if logornot
        set(cb{1},'yticklabel',round(exp(linspace(minV,maxV,cbticks))))
        %         set(cb{2},'yticklabel',round(exp(get(cb{2},'ytick'))))
    end%
    
    %%
%     titalt=['defl-' colorfield '-alt'];
%     savefig(DD.path.plots,200,1200,800,titalt,'dpdf',DD2info(DD));
    %%
    dims=[8*12 6*12];FS=32;
    set(gcf,'paperunits','inch','papersize',dims,'paperposition',[0 0 dims]);
    set(findall(gcf,'type','text'),'FontSize',FS,'interpreter','latex')
    set(gca,'FontSize',FS);
    fname = [DD.path.plots tit];
    print('-depsc','-painters',fname)
    print('-depsc','-painters',fname)
    system(['epstopdf ' fname '.eps']);
    system(['pdfcrop --margins "1 1 1 1" ' fname '.pdf ' fname '.pdf']);
    cpPdfTotexMT(tit);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minV, maxV]=drawColorLinez(ticks,files,fieldName,minlen,cticks,logornot,globe,fac,cmap)
    persistent pp
    if nargin<8,fac=1;end
    %     cmap=jet;% Generate range of color indices that map to cmap
    %% get extremata
    minV=ticks.(cticks)(1);
    maxV=ticks.(cticks)(2);
    if logornot
        minV = log(minV);
        maxV = log(maxV);
    end
    kk=linspace(minV,maxV,size(cmap,1));
    Tac=disp_progress('init','blubb');
    for ee=1:1:numel(files)
        Tac=disp_progress('calc',Tac,round(numel(files)),1000);
        pp=strfind(files{ee},'-d');
        lll = str2double(files{ee}(pp+2:pp+5));
        if lll < minlen,continue,end
        V=load(files{ee},fieldName,'lat','lon');
        VV=V.(fieldName)*fac;
        if logornot,VV = log(VV);end
        cm = spline(kk,cmap',VV);       % Find interpolated colorvalue
        cm(cm>1)=1; cm(cm<0)=0; cm(:,1)=cm(:,2);% Sometimes iterpolation gives values that are out of [0,1] range...
        %% deg2km
        if globe
            yy=V.lat; xx=wrapTo180(V.lon); dJump=100;
        else
            yy=[0 cumsum(deg2km(diff(V.lat)))];
            xx=[0 cumsum(deg2km(diff(V.lon)).*cosd((V.lat(1:end-1)+V.lat(2:end))/2))];
            dJump=1000;
        end
        for ii=1:length(xx)-1
            if  abs(xx(ii+1)-xx(ii))<dJump % avoid 0->360 jumps
                line([xx(ii) xx(ii+1)],[yy(ii) yy(ii+1)],'color',cm(:,ii),'linewidth',0.1,'clipping','off');
            end
        end
    end
    %     caxis([minV maxV])
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUT]=inits(DD)
    disp(['loading maps'])
    OUT.maps=load([DD.path.analyzed.name, 'maps.mat']);
    OUT.la=OUT.maps.Cycs.lat;
    OUT.lo=OUT.maps.Cycs.lon;
    if DD.switchs.netUstuff
        OUT.maps.meanU=load(DD.path.meanU.file);
    end
    %% collect tracks
    OUT.tracksfile=[DD.path.analyzed.name , 'tracks.mat' ];
    
    for ss=1:2
        sen = DD.FieldKeys.senses{ss};
        root=DD.path.analyzedTracks.(sen).name;
        OUT.(sen)={DD.path.analyzedTracks.(sen).files.name};
        Tac=disp_progress('init',['collecting all ', sen]);
        for ff=1:numel(OUT.(sen))
            Tac=disp_progress('calc',Tac,numel(OUT.(sen)),50);
            OUT.tracks.(sen)(ff)={[root OUT.(sen){ff}]};
        end
    end
    %%
    
    %% get vectors
    %     disp(['loading vectors'])
    % 	OUT.vecs=load([DD.path.analyzed.name, 'vecs.mat']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%