%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Sep-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S09_plotsNew
%         DD=initialise([],mfilename);
%         save DD
    load DD
    ticks.rez=get(0,'ScreenPixelsPerInch');
    ticks.width=2*600;
    ticks.height=2*200;
    geo=DD.map.window.geo;
    %     ticks.y= round(linspace(geo.south,geo.north,5));
    ticks.y= [-70 -50 -30 0 30 50 70];
    %     ticks.x=  round(linspace(geo.west,geo.east,5));
    ticks.x=  round(linspace(-180,180,5));
    ticks.axis=[geo.west  geo.east geo.south geo.north];
    ticks.age=[1,5*365,10];
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
    ticks.lat=[ticks.axis(3:4),5];
    %     ticks.minMax=cell2mat(extractfield( load([DD.path.analyzed.name, 'vecs.mat']), 'minMax'));
    %%
    main(DD,ticks)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,T)
   %%
    II=initStuff(DD);
    save S09main II DD T
    %%
    
    sub09_trackstuff
%     %%
%     sub09_mapStuff
%     %%
%     [procData]=inits(DD);
%     save([DD.path.analyzed.name 'procData.mat'],'procData');
%     %     load([DD.path.analyzed.name 'procData.mat'],'procData');
%     %
%     senses = DD.FieldKeys.senses';
%     %   TODO TPz(DD,ticks,procData.tracks,sen,'lat',30,'lat',0);
%     %   TODO TPz(DD,ticks,procData.tracks,sen,'peakampto_mean',30,'amp',1);
%     TPzGlobe(DD,ticks,procData.tracks,senses,'age',50,'age',1,1);
%     %%
%     sub09_histStuff(DD) 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function II=initStuff(DD)
    II.aut=autumn(100);
    II.win=winter(100);
    II.maps=load([DD.path.analyzed.name, 'maps.mat']);  % see S06
    II.la=II.maps.Cycs.lat;
    II.lo=II.maps.Cycs.lon;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=sub09_histStuff(DD)
    
    T(2).files=struct;
    T(1).files = dir2([DD.path.tracks.name '*' DD.FieldKeys.senses{1} '*']);
    T(2).files = dir2([DD.path.tracks.name '*' DD.FieldKeys.senses{2} '*']);
    
    for s=1:2
        N = numel(T(s).files);
        T(s).age=nan(1,N);
        for n=1:N
            T(s).age(n) = str2double(T(s).files(n).name(35:38));
        end
        
        uni = unique(T(s).age);
        H{s} = histc(T(s).age,uni);
    end
    %%
    bar(uni,reshape(cell2mat(H),numel(H{1}),[]),'stacked')
    legend(['anti-cyclones (' num2str(sum(H{1})) ' total)'],['cyclones (' num2str(sum(H{2})) ' total)'])
    axis tight
    xlabel('age [d]')
    ylabel('count')
    set(gca,'xtick',round(linspace(min(uni),max(uni),5)))
    set(gca,'ytick',round(linspace(min(sum(cat(1,H{:}))),max(sum(cat(1,H{:}))),5)))
    %%
    savefig(DD.path.plots,200,800,500,['histTrackCount'],'dpdf',DD2info(DD));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function TPzGlobe(DD,ticks,tracks,senses,colorfield,minlen,cticks,logornot,fac)
    close all
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
    for ss=1:2
        set(cb{ss},'ytick',linspace(ticks.(cticks)(1),ticks.(cticks)(2),ticks.(cticks)(3)))
        if logornot
            set(cb{ss},'yticklabel',round(exp(get(cb{ss},'ytick'))))
        end
    end
    set(cb{2},'yticklabel',[])
    grid minor
    savefig(DD.path.plots,ticks.rez,1400,600,tit,'dpdf')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TPz(DD,ticks,tracks,sen,colorfield,minlen,cticks,logornot)
    drawColorLinez(ticks,tracks.(sen),colorfield,minlen,cticks,logornot,0) ;
    axis([-5000 3000 -2000 2000])
    axis equal
    set(gca,'ytick',[-1000   0   1000])
    set(gca,'xtick',[-4000 -2000 -1000  0 1000])
    tit=['defl-' colorfield '-' sen];
    cb=colorbar;
    if logornot
        set(cb,'yticklabel',round(exp(get(cb,'ytick'))))
    end
    axis tight
    saveas(gcf,[DD.path.plots tit])
    savefig(DD.path.plots,100,3*800,3*500,tit,'dpdf')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
