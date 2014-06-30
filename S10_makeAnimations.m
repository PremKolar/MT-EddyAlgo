%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Apr-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S10_makeAnimations
    DD=initialise;
    %%	set ticks here!
    %     ticks.rez=200;
    ticks.rez=get(0,'ScreenPixelsPerInch');
    %           ticks.rez=42;
    ticks.width=297/25.4*ticks.rez*1;
    ticks.height=ticks.width * DD.map.out.Y/DD.map.out.X;
    %         ticks.height=ticks.width/sqrt(2); % Din a4
    ticks.y= 0;
    ticks.x= 0;
    ticks.age=[1,3*365,10];
    %     ticks.isoper=[DD.thresh.shape.iq,1,10];
    ticks.isoper=[.6,1,10];
    ticks.radius=[20,150,9];
    ticks.radiusToRo=[0.2,5,11];
    ticks.amp=[1,20,7];
    %ticks.visits=[0,max([maps.AntiCycs.visitsSingleEddy(:); maps.Cycs.visitsSingleEddy(:)]),5];
    ticks.visits=[1,20,11];
    ticks.visitsunique=[1,3,3];
    ticks.dist=[-1500;1000;21];
    %ticks.dist=[-100;50;16];
    ticks.disttot=[1;2000;10];
    ticks.vel=[-30;20;6];
    ticks.axis=[DD.map.out.west DD.map.out.east DD.map.out.south DD.map.out.north];
    ticks.lat=[ticks.axis(3:4),5];
    ticks.minMax=cell2mat(extractfield( load([DD.path.analyzed.name, 'vecs.mat']), 'minMax'));
    
    animas(DD)
    %%
    conclude(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function animas(DD)
    file1=[DD.path.eddies.name DD.path.eddies.files(1).name];
    grid=load(cell2mat(extractdeepfield(load(file1),'filename.cut')));
    d.LON=grid.grids.lon;
    d.LAT=grid.grids.lat;
    d.lon=downsize(d.LON,DD.map.out.X*2,DD.map.out.Y);
    d.lat=downsize(d.LAT,DD.map.out.X*2,DD.map.out.Y);
    d.climssh.min=nanmin(grid.grids.ssh(:));
    d.climssh.max=nanmax(grid.grids.ssh(:));
    d.p=[DD.path.plots 'mpngs/'];
    mkdirp(d.p);
    
    parfor ee=1:numel(DD.path.eddies.files)
        disp(num2str(100*ee/numel(DD.path.eddies.files)))
        savepng4mov(d,ee,DD)
    end
    pn=pwd;
    cd(d.p)
    system(['mencoder "mf://flat*.png" -mf fps=10 -o flat.avi -ovc lavc -lavcopts vcodec=mpeg4'])
    system(['mencoder "mf://surf*.png" -mf fps=10 -o surf.avi -ovc lavc -lavcopts vcodec=mpeg4'])
    system(['mplayer  surf.avi'])
    system(['mplayer  flat.avi'])
    cd(pn)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savepng4mov(d,ee,DD)
    d.file=[DD.path.eddies.name DD.path.eddies.files(ee).name];
    ed=load(d.file);
    coor.c=[ed.cyclones.coordinates];
    coor.ac=[ed.anticyclones.coordinates];
    z.c=[ed.cyclones.level];
    z.ac=[ed.anticyclones.level];
    grid=load(cell2mat(extractdeepfield(load(d.file),'filename.cut')));
    ssh=downsize(grid.grids.ssh,2*DD.map.out.X,DD.map.out.Y);
    %%
    pcolor(d.lon,d.lat,ssh);
    caxis([d.climssh.min d.climssh.max]);
    hold on
    sen={'ac';'c'};
    col=[1 0 0; 1 1 1];
    for ss=1:2
        s=sen{ss};
        for cc=1:numel(coor.(s))
            linx=drop_2d_to_1d(coor.(s)(cc).int.y,coor.(s)(cc).int.x,size(d.LAT,1));
            lo=d.LON(linx);
            la=d.LAT(linx);
            plot(lo,la,'color',col(ss,:),'linewidth',1.5);
        end
    end
    
    dayNow.n=ed.anticyclones(1).daynum;
    dayNow.s=datestr(dayNow.n,'yyyymmdd');
    for tt=1:numel(DD.path.tracks.files)
        from=datenum(DD.path.tracks.files(tt).name(6:13),'yyyymmdd');
        till=datenum(DD.path.tracks.files(tt).name(15:22),'yyyymmdd');
        if from<=dayNow.n && till>=dayNow.n
            t.file=[DD.path.tracks.name DD.path.tracks.files(tt).name];
            tr=load(t.file); tr=tr.trck;
            daynums=cat(2,tr.daynum);
            [~,nn]=min(abs(daynums-dayNow.n));
            tr=tr(1:nn);
            lo=extractdeepfield(tr,'geo.lon');
            la=extractdeepfield(tr,'geo.lat');
            plot(lo,la,'color',rainbow(1,1,1,max([100-nn,0]),100),'linewidth',1.5);
        end
    end
    
    
    savefig2png4mov(d.p,100,1600,1200,sprintf('flat%06d',ee))
    
    %%
    figure
    surf(d.lon,d.lat,ssh);
    axis([min(d.lon(:)) max(d.lon(:)) min(d.lat(:)) max(d.lat(:)) d.climssh.min*2 d.climssh.max*2 d.climssh.min/2 d.climssh.max/2]);
    hold on
    sen={'ac';'c'};
    col=[1 0 0; 1 1 1];
    for ss=1:2
        s=sen{ss};
        for cc=1:numel(coor.(s))
            linx=drop_2d_to_1d(coor.(s)(cc).int.y,coor.(s)(cc).int.x,size(d.LAT,1));
            lo=d.LON(linx);
            la=d.LAT(linx);
            zz=repmat(z.(s)(cc),1,numel(lo));
            plot3(lo,la,zz,'color',col(ss,:),'linewidth',2);
        end
    end
    
    
    
    
    
    savefig2png4mov(d.p,100,1600,1200,sprintf('surf%06d',ee))
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
