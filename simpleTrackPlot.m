function simpleTrackPlot
    DD = initialise([],mfilename);
    save DD
%     system(['cp DD.mat ../aviSingleTrackI/DDII.mat'])
    addpath(genpath('./'))
    %     load DD
    senses = DD.FieldKeys.senses;
    ss=2;
    SSs=[-1 1];
    sen = senses{ss};
    %%
    xx = [20 100];
    yy = [20 60];
    %%
    kk = 0;
    tracks = dir2([DD.path.tracks.name]);
    tracks(1:2)  = [];
    TRACKS = loadtracks(tracks);
%     system(['cp TRACKS.mat ../aviSingleTrackI/TRACKSII.mat'])
    
        %%
        TT = DD.time.from.num:DD.time.delta_t:DD.time.till.num;
           for  tt = TT
            kk=kk+1;
            main(tt,kk,DD,xx,yy,SSs,ss,sen,TRACKS);
        end
    
end

function main(tt,kk,DD,xx,yy,SSs,ss,sen,TRACKS)
    clf
    clc;
    fprintf('day %d/%d\n',kk,ceil((DD.time.till.num-DD.time.from.num)/DD.time.delta_t));
    
    cut = load([DD.path.cuts.name DD.path.cuts.files(kk).name]);
    eddy =load([DD.path.eddies.name DD.path.eddies.files(kk).name]);
    ssh = cut.fields.ssh(yy(1):yy(2),xx(1):xx(2));
    lon = cut.fields.lon(yy(1):yy(2),xx(1):xx(2));
    lat = cut.fields.lat(yy(1):yy(2),xx(1):xx(2));
    %%
    plotssh(lon,lat,ssh)
    %%
    ploteddies(eddy,lon,lat,sen,xx,yy)
    %%
    plottracks(TRACKS,tt,lon,lat,SSs(ss),xx,yy);
    %%
    %     title(datestr(tt),'fontsize',18,'interpreter','latex');
    set(gca,'xtick', 40:.25: 80,'fontsize',2)
    set(gca,'ytick',-50:.25:-30,'fontsize',2)
    %     set(gca,'
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')
    set(gca,'layer','top')
    grid on
    
    %%
    saveas(gcf,sprintf('%03d.png',kk))
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotssh(lon,lat,ssh)
    pcolor(lon,lat,ssh);
    shading flat
    shading interp
    hold on;
    %     colormap([bone(30);flipud(hot(30))]);
    colormap(parula(50));
    caxis([-.7,.4])
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploteddies(eddy,lon,lat,sen,xx,yy)
    for ee = 1:numel(eddy.(sen))
        coor = eddy.(sen)(ee).coor.exact;
        x = coor.x - xx(1) + 1;
        y = coor.y - yy(1) + 1;
        lo = interp2(lon,x,y);
        la = interp2(lat,x,y);
        plot(lo,la,'black','linestyle','--')
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plottracks(T,tt,lon,lat,senN,xx,yy)
    
    for jj = 1:numel(T)
        %         fprintf('plotting track %d/%d\n',jj,numel(T));
        dateA =  T(jj).d.track(1).daynum;
        dateB =  T(jj).d.track(end).daynum;
        
        if dateA <= tt && dateB >= tt && T(jj).d.track(1).sense.num == senN
            track = T(jj).d.track;
            tdn = cat(1,track.daynum);
            toplot = tdn <= tt   ;
            x = extractdeepfield(track(toplot),'trackref.x')- xx(1) + 1;
            y = extractdeepfield(track(toplot),'trackref.y')- yy(1) + 1;
            lo = interp2(lon,x,y);
            la = interp2(lat,x,y);
            plot(lo,la,'black','linewidth',2);
        end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TRACKS = loadtracks(tracks)
    if ~exist('TRACKS.mat','file')
        TRACKS(numel(tracks)) = struct;
        for jj = 1:numel(tracks)
            clc;
            fprintf('loading %d/%d\n',jj,numel(tracks));
            TRACKS(jj).d = load(tracks(jj).fullname);
        end
        save TRACKS TRACKS
    else
        load TRACKS
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%