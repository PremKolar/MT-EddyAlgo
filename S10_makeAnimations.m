%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Apr-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S10_makeAnimations
    DD=initialise([],mfilename);
    ticks.rez=get(0,'ScreenPixelsPerInch');
    ticks.width=297/25.4*ticks.rez*1;
    ticks.height=ticks.width * DD.map.out.Y/DD.map.out.X;
    ticks.y= 0;
    ticks.x= 0;
    ticks.age=[1,3*365,10];
    ticks.isoper=[.6,1,10];
    ticks.radius=[20,150,9];
    ticks.radiusToRo=[0.2,5,11];
    ticks.amp=[1,20,7];
    ticks.visits=[1,20,11];
    ticks.visitsunique=[1,3,3];
    ticks.dist=[-1500;1000;21];
    ticks.disttot=[1;2000;10];
    ticks.vel=[-30;20;6];
    ticks.axis=[DD.map.out.west DD.map.out.east DD.map.out.south DD.map.out.north];
    ticks.lat=[ticks.axis(3:4),5];
    ticks.minMax=cell2mat(extractfield( load([DD.path.analyzed.name, 'vecs.mat']), 'minMax'));
    animas(DD)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function animas(DD)
    file1=[DD.path.eddies.name DD.path.eddies.files(1).name];
    
    grid=load(cell2mat(extractdeepfield(load(file1),'filename.cut')));
    d.LON=grid.grids.lon;
    d.LAT=grid.grids.lat;
    d.lon=grid.grids.lon;
    d.lat=grid.grids.lat;
    d.climssh.min=nanmin(grid.grids.ssh(:));
    d.climssh.max=nanmax(grid.grids.ssh(:));
    d.p=[DD.path.plots 'mpngs/'];
    mkdirp(d.p);
    frms=3800;
    range=round(linspace(1,numel(DD.path.eddies.files),frms));
    for cc=1:numel(range)
        ee=range(cc)%
        disp(num2str(100*ee/numel(DD.path.eddies.files)));
        savepng4mov(d,ee,DD)
    end
    fps=min([1 round(frms/60)]);
    pn=pwd;
    cd(d.p)
    system(['mencoder "mf://flat*.png" -mf fps=' num2str(fps) ' -o flat.avi -ovc lavc -lavcopts vcodec=mpeg4'])
    system(['mplayer  flat.avi'])
    cd(pn)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savepng4mov(d,ee,DD)
    d.file=[DD.path.eddies.name DD.path.eddies.files(ee).name];
    if exist([ d.p sprintf('flat%06d.png',ee)],'file')
        return
    end
    [~,fn,~] = fileparts(d.file);
    d.dtnm=datenum(fn(7:14),'yyyymmdd');
    ed=load(d.file);
    coor.c=[ed.cyclones.coordinates];
    coor.ac=[ed.anticyclones.coordinates];
    z.c=[ed.cyclones.level];
    z.ac=[ed.anticyclones.level];
    try
        grid=load(cell2mat(extractdeepfield(load(d.file),'filename.cut')));
    catch
        return
    end
    ssh=grid.grids.ssh;
    pcolor(d.lon,d.lat,ssh);
    shading flat
    axis equal tight
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
    
    title([num2mstr(ee)])
    xlabel(datestr(d.dtnm))
    savefig2png4mov(d.p,100,800,600,sprintf('flat%06d',ee))
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
