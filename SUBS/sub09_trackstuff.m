function sub09_trackstuff
    load S09main II DD T
%     sub09_trackinit(DD);
    TR=getTR(DD) ;
    %%
    senses=DD.FieldKeys.senses;
    catsen= @(f) [TR.(senses{1}).(f); TR.(senses{2}).(f) ];
    S.t2l=@(t) round(linspace(t(1),t(2),t(3)));
    %%
    rad=round(catsen('rad')/1000);
    radLe=round(catsen('radLe')/1000);
    radLeff=round(catsen('radLeff')/1000);
    radL=round(catsen('radL')/1000);
    vel=catsen('vel')*100;
    age=catsen('age');
    lat=catsen('lat');
    lon=catsen('lon'); %#ok<NASGU>
   
%     reflin=catsen('reflin');
    %%
    S.rightyscalenum=5;
    age(end+1:end+S.rightyscalenum)=max(age)-0;
    lat(end+1:end+S.rightyscalenum)=S.t2l([min(lat) max(lat) S.rightyscalenum]);
    lon(end+1:end+S.rightyscalenum)=S.t2l([min(lon) max(lon) S.rightyscalenum]);
    rad(end+1:end+S.rightyscalenum)=S.t2l([min(rad) max(rad) S.rightyscalenum]);
    vel(end+1:end+S.rightyscalenum)=10;

    radL(end+1:end+S.rightyscalenum)=0;
    radLe(end+1:end+S.rightyscalenum)=0;
    radLeff(end+1:end+S.rightyscalenum)=0;
%     reflin(end+1:end+S.rightyscalenum)=0;

    %%
    [~,sml2lrg] = sort(rad)  ;
    S.age = age(fliplr(sml2lrg));
    S.lat = lat(fliplr(sml2lrg));
    S.lon = lon(fliplr(sml2lrg));
    S.rad = rad(fliplr(sml2lrg));
    S.vel = vel(fliplr(sml2lrg));

    S.radLe   = radLe(fliplr(sml2lrg));
    S.radLeff = radLeff(fliplr(sml2lrg));
    S.radL    = radL(fliplr(sml2lrg));

%     S.reflin  = reflin(fliplr(sml2lrg));

    %% kill unrealistic data
    zerage  = S.age<=0  ;
    velHigh = S.vel>20 | S.vel <-30;
    radnill = isnan(S.rad) | S.rad==0;
    killTag = zerage | velHigh | radnill ;

    FN=fieldnames(S);
    for ii=1:numel(FN)
        try
            S.(FN{ii})(killTag)=[];
        end
    end

    %%
    spmdblock(S,DD,II,T);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmdblock(S,DD,II,T)
    velZonmeans(S,DD,II,T);
    scaleZonmeans(S,DD,II,T);
    %             	scattStuff(S,T,DD,II);

    %     	spmd
    %     		switch labindex
    %     			case 1
    %     				scaleZonmeans(S,DD,II,T);
    %     			case 2
    %     				velZonmeans(S,DD,II,T);
    % %     			case 3
    % %     				scattStuff(S,T,DD,II);
    %     		end
    %      	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=scaleZonmeans(S,DD,II,T) %#ok<INUSD>
    close all
    chelt = imread('/scratch/uni/ifmto/u300065/FINAL/presStuff/LTpresMT/FIGS/png1024x/chSc.png');
    LA     = round(S.lat);
    LAuniq = unique(LA)';
    %     FN     = {'rad','radL','radLe','radLeff'};
    FN     = {'rad','radLe','radLeff'};
    %     FN     = {'Lrossby'};

%     Rpath = DD.path.Rossby.name;
%     Rname = [DD.FieldKeys.Rossby{2} ,'.mat'];
%     LR = getfield(load([Rpath Rname]),'data');
%     zerFlag = S.reflin == 0;
%     S.reflin(zerFlag) = 1;
%     S.Lrossby = LR(S.reflin)/1000; % m2km
%     S.Lrossby(zerFlag) = nan;
    for ff=1:numel(FN)
        fn=FN{ff}
        S.(fn)(S.(fn)<10) = nan;
    end
    %%
    visits = nan(size(LAuniq));
    for ff=1:numel(FN)
        fn=FN{ff}
        vvM(numel(LAuniq)).(fn)=struct;
        [vvM(:).(fn)]=deal(nan);

        for cc=1:(numel(LAuniq))
            idx=LA==LAuniq(cc);
            visits(cc) = sum(idx);

            if visits(cc) >= 100
                vvM(cc).(fn)=nanmedian(S.(fn)(idx));
                if abs(LAuniq(cc))<=5
                    vvM(cc).(fn)=nan;
                end

            end
        end
    end



    %%
    h.ch=chOverLayScale(chelt,LAuniq,vvM);
    savefig(DD.path.plots,100,800,800,['S-scaleZonmean4chelt11comp'],'dpdf',DD2info(DD));


    %     legend('off')
    %     title([''])
    %     savefig(DD.path.plots,100,800,800,['S-RossbyLfromPopToCh'],'dpdf',DD2info(DD));
    %%
    % 	[h.own,pp,dd]=ownPlotScale(DD,II,LAuniq,vvM,vvS); %#ok<NASGU>
    % 	[~,pw]=fileparts(pwd);
    % 	save(sprintf('scaleZonMean-%s.mat',pw),'h','pp','dd');
    % 	savefig(DD.path.plots,100,800,800,['S-scaleZonmean'],'dpdf',DD2info(DD));

    %% TODO
    figure(2)
    cc = 70;
    fn=FN{1}
    idx=LA==LAuniq(cc); % -10
    hist(S.(fn)(idx),50)
    axis tight
    xlabel('\sigma at -10^{\circ}')
    title(sprintf('total: %d counts',sum(idx)))
    savefig(DD.path.plots,100,600,600,['hist-sigmaAt-10deg'],'dpdf',DD2info(DD));


    figure(3)
    cc = 30;
    fn=FN{1}
    idx=LA==LAuniq(cc); % -50
    hist(S.(fn)(idx),50)
    axis tight
    title(sprintf('total: %d counts',sum(idx)))
    xlabel('\sigma at -50^{\circ}')
    savefig(DD.path.plots,100,600,600,['hist-sigmaAt-50deg'],'dpdf',DD2info(DD));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=velZonmeans(S,DD,II,T) %#ok<INUSD>
    close all

%     Rpath = DD.path.Rossby.name;
%     Rname = [DD.FieldKeys.Rossby{1} ,'.mat'];
%     cR = getfield(load([Rpath Rname]),'data');
%     zerFlag = S.reflin == 0;
%     S.reflin(zerFlag) = 1;
%     S.Crossby = cR(S.reflin)*100; % m2cm
%     S.Crossby(zerFlag) = nan;


  LA     = round(S.lon);
    LAuniq = unique(LA)';
    vvM=nan(size(LAuniq));
   
   
    visits = nan(size(LAuniq));
    for cc=1:(numel(LAuniq))
        idx=LA==LAuniq(cc);
        visits(cc) = sum(idx);
    end








    close all
    LA     = round(S.lat);
    LAuniq = unique(LA)';
    vvM=nan(size(LAuniq));
    vvS=nan(size(LAuniq));

    vvSkew=nan(size(LAuniq));
    visits = nan(size(LAuniq));
    for cc=1:(numel(LAuniq))
        idx=LA==LAuniq(cc);
        visits(cc) = sum(idx);
        %         if visits(cc) >= 100
        vvM(cc)=nanmedian(S.vel(idx));
        vvS(cc)=std(S.vel(idx));
        vvSkew(cc)=skewness(S.vel(idx));
        %         end
    end
    %%
    % TODO do this with pop7 or better pop3 data ! and maybe do similar with
    % scales..
    figure(10)
    clf
    hold on
    plot(LAuniq,smooth(-vvM,10),'blue')
    %  plot(LAuniq,vvS,'red')
    plot(LAuniq,smooth(-vvSkew*5,10),'green')
    plot(LAuniq,smooth(visits/20000*10,10),'red')
    axis([-80 80 -2 10])
    legend({'-u','-skewness(u)','count'})
    plot([-80 80],[0 0],'--black')
    grid on
    set(gca,'yticklabel','')
    savefig(DD.path.plots,100,800,200,['Skew'],'dpdf',DD2info(DD));


    vvM(abs(LAuniq)<5)=nan;
    vvS(abs(LAuniq)<5)=nan;


    %%
    %     [h.own,~,dd]=ownPlotVel(DD,II,LAuniq,vvM,vvS); %#ok<NASGU>
    %     [~,pw]=fileparts(pwd);
    %     save(sprintf('velZonMean-%s.mat',pw),'h','pp','dd');
    %     savefig(DD.path.plots,100,800,800,['S-velZonmean'],'dpdf',DD2info(DD));
    %%
    chelt = imread('/scratch/uni/ifmto/u300065/FINAL/PLOTS/chelt11Ucomp.jpg');
    chelt= chelt(135:3595,415:3790,:);
    h.ch=chOverLay(S,DD,chelt,LAuniq,vvM);
    savefig(DD.path.plots,100,800,800,['S-velZonmean4chelt11comp'],'dpdf',DD2info(DD));

    %       figure
    %     h.ch=chOverLay(S,DD,chelt,LAuniq,vvCross);
    %     title([])
    %     savefig(DD.path.plots,100,800,800,['S-RossbyCfromPopToCh'],'dpdf',DD2info(DD));



    %% TODO
    figure(2)
    cc = 70;
    idx=LA==LAuniq(cc); % -10
    hist(S.vel(idx),50)
    axis tight
    xlabel('u [cm/s] at -10^{\circ}')
    title(sprintf('total: %d counts',sum(idx)))
    savefig(DD.path.plots,100,600,600,['hist-uAt-10deg'],'dpdf',DD2info(DD));


    figure(3)
    cc = 30;
    idx=LA==LAuniq(cc); % -50
    hist(S.vel(idx),50)
    axis tight
    xlabel('u [cm/s] at -50^{\circ}')
    title(sprintf('total: %d counts',sum(idx)))
    savefig(DD.path.plots,100,600,600,['hist-uAt-50deg'],'dpdf',DD2info(DD));



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,pp,dd]=ownPlotVel(DD,II,LAuniq,vvM,vvS)
    %%
    clf
    lw=2;
    %%
    dd(1).y=-II.maps.zonMean.Rossby.small.phaseSpeed*100;
    dd(2).y=-vvM;
    dd(4).y=-vvM+vvS;
    dd(5).y=-vvM-vvS;
    dd(3).y= [0 0];
    %%
    dd(1).name='rossbyPhaseSpeed';
    dd(2).name='zonal mean eddy phase speeds';
    dd(4).name='std upper bound';
    dd(5).name='std lower bound';
    dd(3).name='nill line';
    %%
    dd(1).x=II.la(:,1);
    dd(2).x=LAuniq;
    dd(4).x=LAuniq;
    dd(5).x=LAuniq;
    geo=DD.map.window.geo;
    dd(3).x=[geo.south geo.north];
    %%
    pp(1)=plot(dd(1).x,dd(1).y); 	hold on
    pp(2)=plot(dd(2).x,dd(2).y,'r');
    pp(4)=plot(dd(4).x,dd(4).y,'r');
    pp(5)=plot(dd(5).x,dd(5).y,'r');
    pp(3)=plot(dd(3).x,dd(3).y,'b--');
    %%
    axis([-70 70 -5 20])
    set(pp(1:3),'linewidth',lw)
    leg=legend('Rossby-wave phase-speed',2,'all eddies',2,'std');
    legch=get(leg,'children');
    set( legch(1:3),'linewidth',lw)
    ylabel('[cm/s]')
    xlabel('[latitude]')
    title(['westward propagation [cm/s]'])
    set(get(gcf,'children'),'linewidth',lw)
    grid on
    h=gcf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,pp,dd]=ownPlotScale(DD,II,LAuniq,vvM,vvS)
    lw=2;
    %%
    dd(1).y=II.maps.zonMean.Rossby.small.radius/1000*2;
    dd(2).y=vvM;
    dd(4).y=vvM+vvS;
    dd(5).y=vvM-vvS;
    dd(3).y= [0 0];
    %%
    dd(1).name='2x 1st barocl. Rossby radius';
    dd(2).name='zonal mean eddy radius (\sigma)';
    dd(4).name='std upper bound';
    dd(5).name='std lower bound';
    dd(3).name='nill line';
    %%
    dd(1).x=II.la(:,1);
    dd(2).x=LAuniq;
    dd(4).x=LAuniq;
    dd(5).x=LAuniq;
    geo=DD.map.window.geo;
    dd(3).x=[geo.south geo.north];
    %%
    clf
    pp(1)=plot(dd(1).x,dd(1).y); 	hold on
    pp(2)=plot(dd(2).x,dd(2).y,'r');
    pp(4)=plot(dd(4).x,dd(4).y,'r');
    pp(5)=plot(dd(5).x,dd(5).y,'r');
    pp(3)=plot(dd(3).x,dd(3).y,'b--');
    %%
    axis([-70 70 0 max(dd(4).y)])
    set(pp(1:3),'linewidth',lw)
    leg=legend(dd(1).name,2,dd(2).name,2,'std');
    legch=get(leg,'children');
    set( legch(1:3),'linewidth',lw)
    ylabel('[km]')
    xlabel('[latitude]')
    set(get(gcf,'children'),'linewidth',lw)
    grid on
    h=gcf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=chOverLayScale(chelt,LAuniq,vvM)
    clf
    [Y,X,Z]=size(chelt);
    ch=reshape(flipud(reshape(chelt,Y,[])),[Y,X,Z]);
    ch=permute(reshape(fliplr(reshape(permute(ch,[3,1,2]),1,[])),[Z,Y,X]),[2,3,1]);
    ch=permute(reshape(flipud(reshape(permute(ch,[2,1,3]),[X,Y*Z])),[X,Y,Z]),[2,1,3]);
    %%
    imagesc(linspace(-70,70,X),linspace(0,275,Y),ch)
    hold on
    %%
    FN=fieldnames(vvM)';

    for ff=1:numel(FN)
        fn=FN{ff};
        lau=LAuniq;
        kill.lau = isnan(lau) | abs(lau)>70;
        vvm=-cat(2,vvM.(fn))+275;
%         vvm=-vvM.(fn)+275;
        kill.vvm = isnan(vvm);
        kill.both = kill.lau | kill.vvm;
        vvm(kill.both)=nan;
        lau(kill.both)=nan;
        mid=floor(numel(lau)/2);
        x.a=lau(1:mid);
        x.b=lau(mid+2:end);
        y.a=vvm(1:mid);
        y.b=vvm(mid+2:end);
        forleg(ff)=plot(x.a,y.a,'color',int2color(ff),'linewidth',1);
        plot(x.b,y.b,'color',int2color(ff),'linewidth',1);
        plot(lau,vvm,'color',int2color(ff),'linestyle','.','markersize',8);
    end
    legend(forleg,FN)
    %%
    set(gca, 'ytick', 0:25:275);
    axis([-70 70 0 275])
    set(gca, 'yticklabel', flipud(get(gca,'yticklabel')));
    axis tight
    grid on
    ylabel('[km]')
    xlabel('[latitude]')
    title(['\sigma'])
    h.fig=gcf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=chOverLay(S,DD,chelt,LAuniq,vvM)
    [Y,X,Z]=size(chelt);
    ch=reshape(flipud(reshape(chelt,Y,[])),[Y,X,Z]);
    ch=permute(reshape(fliplr(reshape(permute(ch,[3,1,2]),1,[])),[Z,Y,X]),[2,3,1]);
    ch=permute(reshape(flipud(reshape(permute(ch,[2,1,3]),[X,Y*Z])),[X,Y,Z]),[2,1,3]);
    %%
    vvm=vvM+15;
    lau=LAuniq;
    kill=isnan(lau) | isnan(vvM) | abs(lau)<10 | abs(lau)>65;
    lau(kill)=nan;
    vvm(kill)=nan;
    %%
    clf
    imagesc(linspace(-50,50,X),linspace(-5,20,Y),ch)
    hold on
    mid=floor(numel(lau)/2);
    x.a=lau(1:mid);
    x.b=lau(mid+2:end);
    y.a=vvm(1:mid);
    y.b=vvm(mid+2:end);
    plot(x.a,y.a,'g-','linewidth',1);
    plot(x.b,y.b,'g-','linewidth',1);
    plot(lau,vvm,'g.','markersize',8);
    set(gca, 'yticklabel', flipud(get(gca,'yticklabel')));
    axis tight
    grid on
    ylabel('[cm/s]')
    xlabel('[latitude]')
    title(['westward propagation [cm/s]'])
    axis([-65 65 -5 20])
    h=gcf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scattStuff(S,T,DD,II)
    age=S.age;
    lat=S.lat;
    vel=S.vel;
    rad=round(S.rad);
    inc=1;
    %%
    clf
    oie=@(inc,x) x(1:inc:end);
    incscatter=@(inc,a,b,c,d) scatter(oie(inc,a),oie(inc,b),oie(inc,c),oie(inc,d));
    incscatter(inc,age,lat,rad,vel);
    grid on
    axis tight
    set(gca,'XAxisLocation','bottom')
    cb=colorbar;
    cb1 = findobj(gcf,'Type','axes','Tag','Colorbar');
    cbIm = findobj(cb1,'Type','image');
    alpha(cbIm,0.5)
    set(cb,'location','north','xtick',(S.t2l(T.vel)),'xlim',T.vel([1 2]))
    doublemap([T.vel(1) 0 T.vel(2)],II.aut,II.win,[.9 1 .9],20)
    h1=gca;
    h1pos = get(h1,'Position'); % store position of first axes
    h2 = axes('Position',h1pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none');
    set(h2, ...
        'ytick',linspace(0,1,S.rightyscalenum),...
        'xtick',[],...
        'yticklabel',(S.t2l([min(rad) max(rad) S.rightyscalenum])))
    set(h1, ...
        'ytick',[-50 -30 -10 0 10 30 50],...
        'xtick',S.t2l(T.age))
    ylabel(h2,'radius [km]')
    ylabel(h1,'lat  [$^{\circ}$]')
    xlabel(h1,'age [d]')
    xlabel(h2,'zon. vel.  [cm/s] - eddies beyond scale dismissed!')
    set(get(gcf,'children'),'clipping','off')
    %%
    %  figure(1)
    savefig(DD.path.plots,T.rez,T.width,T.height,['sct-ageLatRadU'],'dpng',DD2info(DD));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TR=getTR(DD)
    xlt=@(sen,f) extractfield(load(['TR-' sen '-' f '.mat']),'tmp');
    %     F={'rad','age','lat','lon'};
    F={'rad','age','lat','lon','radL','radLeff','radLe'};
%     F={'rad','age','lat','lon','radL','radLeff','radLe','reflin'};
    g=@(c) cat(1,c{:});
    for ss=1:2
        for fi=1:numel(F);f=F{fi};
            sen=DD.FieldKeys.senses{ss};
            TR.(sen).(f)=((xlt(sen,f))');
        end
        f='vel';
        TR.(sen).(f)=g(g(xlt(sen,f)));
    end
end
