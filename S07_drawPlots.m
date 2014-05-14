%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Apr-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all figs are saved to ~/FIGS/  !
function S07_drawPlots
%% init

 init_threads(3);
 spmd
    [DD,threadData]=inits;
    
    %%	set ticks here!
    ticks.rez=100;
    ticks.width=297/25.4*ticks.rez;
    ticks.height=ticks.width/sqrt(2);
    % 	ticks.y= linspace(20,60,2);
    % 	ticks.x= linspace(-90,-50,2);
    ticks.y= 0;
    ticks.x= 0;
    ticks.age=[1,DD.time.span/2,10];
    ticks.isoper=[DD.thresh.shape.iq,1,10];
    ticks.radius=[50,150,6];
    ticks.amp=[1,20,7];
    %ticks.visits=[0,max([maps.AntiCycs.visitsSingleEddy(:); maps.Cycs.visitsSingleEddy(:)]),5];
    ticks.visits=[1,20,11];
    ticks.dist=[-2000;1000;7];
    %ticks.dist=[-100;50;16];
    ticks.disttot=[1;2000;14];
    ticks.vel=[-30;20;6];
    ticks.axis=[DD.dim.west DD.dim.east DD.dim.south DD.dim.north];
    ticks.lat=[ticks.axis(3:4),5];
    
    
    %% main
    main(DD,threadData,ticks);
    
 end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,OUT]=inits

DD=initialise();
if labindex==2
    OUT.maps=load([DD.path.analyzed.name, 'maps.mat']);
    OUT.la=OUT.maps.Cycs.GLA;
    OUT.lo=OUT.maps.Cycs.GLO;
end
%% collect tracks
if labindex==1
    OUT.tracksfile=[DD.path.analyzed.name , 'tracks.mat' ];  
        root=DD.path.analyzedTracks.AC.name;
        OUT.ACs={DD.path.analyzedTracks.AC.files.name};
        Tac=disp_progress('init','collecting all ACs');
        for ff=1:numel(OUT.ACs)
            Tac=disp_progress('calc',Tac,numel(OUT.ACs),100);
            OUT.tracks.AntiCycs(ff)={[root OUT.ACs{ff}]};
        end
        %%
        root=DD.path.analyzedTracks.C.name;
        OUT.Cs={DD.path.analyzedTracks.C.files.name};
        Tc=disp_progress('init','collecting all Cs');
        for ff=1:numel(OUT.Cs)
            Tc=disp_progress('calc',Tc,numel(OUT.Cs),100);
            OUT.tracks.Cycs(ff)={[root OUT.Cs{ff}]};
        end     
end

%% get vectors
if labindex==3 || labindex==2
    OUT.vecs=load([DD.path.analyzed.name, 'vecs.mat']);
end
end

function main(DD,IN,ticks)
if labindex==3
    histstuff(IN.vecs,DD,ticks)
end
if labindex==1
    try
    trackPlots(DD,ticks,IN.tracks.tracks)
    catch
        trackPlots(DD,ticks,IN.tracks)
    end
end
if labindex==2
    mapstuff(IN.maps,IN.vecs,DD,ticks,IN.lo,IN.la)
end
end
function histstuff(vecs,DD,ticks)
for sense=fieldnames(vecs)';sen=sense{1};
    figure
    lat=vecs.(sen).lat;
    range.(sen).lat=round(nanmin(vecs.(sen).lat)):round(nanmax(vecs.(sen).lat));
    hc.(sen).lat=histc(lat,range.(sen).lat);
    semilogy(range.(sen).lat,hc.(sen).lat)
    tit=['number of ',sen,' per 1° lat']
    title([tit])
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_latNum'],0);
    %%
    figure
    age=vecs.(sen).age;
    range.(sen).age=round(nanmin(vecs.(sen).age)):round(nanmax(vecs.(sen).age));
    hc.(sen).age=histc(age,range.(sen).age);
    semilogy(range.(sen).age,hc.(sen).age)
    tit=['number of ',sen,' per age'];
    unit='d';
    title([tit,' [',unit,']'])
    hc.(sen).cum.age=fliplr(cumsum(fliplr(hc.(sen).age)));
    semilogy(range.(sen).age,hc.(sen).cum.age)
    tit=['upper tail cumulative of ',sen,' per age'];
    title([tit,' [',unit,']'])
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_ageUpCum'],0);
end
%%
figure
len=nanmin([length(hc.Cycs.lat),length(hc.AntiCycs.lat)]);
hc.rat.lat=hc.Cycs.lat(1:len)./hc.AntiCycs.lat(1:len);
hc.rat.lat(isinf(hc.rat.lat))=nan;
hc.rat.lat((hc.rat.lat==0))=nan;
totnum=hc.Cycs.lat(1:len) + hc.AntiCycs.lat(1:len);
x=(range.Cycs.lat(1:len));
y=log10(hc.rat.lat);
x(totnum==0)=[];
y(totnum==0)=[];
totnum(totnum==0)=[];
colors = flipud(bone(max(totnum)));
cols=colors(totnum,:);
for ii = 1:length(x)
    a=bar(x(ii), y(ii));
    hold on
    set(a,'facecolor', cols(ii,:));
end

axis([range.Cycs.lat(1) range.Cycs.lat(end) -1 1]);
ticks.y=linspace(-1,1,3);
ticks.lab.y=num2cell(10.^(ticks.y));
set(gca,'ytick',ticks.y)
set(gca,'yticklabel',ticks.lab.y)
tit=['ratio of cyclones to anticyclones as function of latitude'];
title(tit)
cm=flipud(bone);
colormap(cm)
cb=colorbar;
zticks=linspace(0,max(totnum),5)';
zticklabel=num2str(round((zticks)));
caxis([zticks(1) zticks(end)])
set(cb,'ytick',zticks);
set(cb,'yticklabel',zticklabel);

savefig(ticks.rez,ticks.width,ticks.height,[sen,'_RatLat'],0);
%%
figure
len=min([length(hc.AntiCycs.age), length(hc.Cycs.age)]);
hc.rat.age=hc.Cycs.age(1:len)./hc.AntiCycs.age(1:len);
hc.rat.age(isinf(hc.rat.age))=nan;
hc.rat.age((hc.rat.age==0))=nan;
totnum=hc.Cycs.age(1:len) + hc.AntiCycs.age(1:len);
x=sqrt(range.Cycs.age(1:len));
y=log10(hc.rat.age);
x(totnum<2)=[];
y(totnum<2)=[];
totnum(totnum<2)=[];
tnlog=log(totnum);
tmin=min(tnlog);
tmax=max(tnlog);
cmap=flipud(bone);
kk=(linspace(tmin,tmax,size(cmap,1)));

for ii = 1:length(x)
    a=bar(x(ii), y(ii),.1);
    hold on
    cm = spline(kk,cmap',(tnlog(ii)));  % Find interpolated colorvalue
    cm(cm<0)=0;
    cm(cm>1)=1;
    set(a,'facecolor', cm');
end

axis([min(x)-1 max(x) -1 1]);
xtck=get(gca,'xtick');
set(gca,'xticklabel',xtck.^2)
ticks.y=linspace(-1,1,3);
ticks.lab.y=num2cell(10.^(ticks.y));
set(gca,'ytick',ticks.y)
set(gca,'yticklabel',ticks.lab.y)
tit='ratio of cyclones to anticyclones as function of age';
xlab='age [d] - gray scale indicates total number of eddies available';
title(tit)
xlabel(xlab)
cm=flipud(bone);
colormap(cm)
cb=colorbar;
zticks=linspace(0,max(totnum),5)';
zticklabel=num2str(round((zticks)));
caxis([zticks(1) zticks(end)])
set(cb,'ytick',zticks);
set(cb,'yticklabel',zticklabel);
savefig(ticks.rez,ticks.width,ticks.height,[sen,'_LatNum'],0);
%%
figure
hc.rat.cum.age=hc.Cycs.cum.age(1:len)./hc.AntiCycs.cum.age(1:len);
hc.rat.cum.age(isinf(hc.rat.cum.age))=nan;
plot(range.Cycs.age(1:len),log10(hc.rat.cum.age));
axis([range.Cycs.age(1) range.Cycs.age(len) -1 1]);
ticks.y=linspace(-1,1,3);
ticks.lab.y=num2cell(10.^(ticks.y));
set(gca,'ytick',ticks.y)
set(gca,'yticklabel',ticks.lab.y)
tit=['ratio of cyclones to anticyclones as upper tail cum. of age'];
title(tit)
savefig(ticks.rez,ticks.width,ticks.height,[sen,'_RatAge'],0);
end

function mapstuff(maps,vecs,DD,ticks,lo,la)
aut=autumn;
win=winter;
for sense=fieldnames(maps)';sen=sense{1};
    
    
    
    figure
    b.la=vecs.(sen).birth.lat;
    b.lo=vecs.(sen).birth.lon;
    d.la=vecs.(sen).death.lat;
    d.lo=vecs.(sen).death.lon;
    plot(b.lo,b.la,'.r',d.lo,d.la,'.g','markersize',5)
    hold on
    drawcoast
    axis(ticks.axis);
    set(gca,'ytick',ticks.y);
    set(gca,'xtick',ticks.x) ;
    legend('births','deaths')
    title(sen)
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_deathsBirths'],0);
    
    %%
    figure
    maps.(sen).age.logmean=log(maps.(sen).age.mean);
    VV=maps.(sen).age.logmean;
    pcolor(lo,la,VV);shading flat
    caxis([ticks.age(1) ticks.age(2)])
    decorate('age',ticks,DD,sen,'age','d',1,1);
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapAge'],0)
    %%
    figure
    VV=maps.(sen).visits.single;
    VV(VV==0)=nan;
    pcolor(lo,la,log(VV));shading flat
    caxis([ticks.visits(1) ticks.visits(2)])
    decorate('visits',ticks,DD,sen,'Visits of unique eddy','',1,1);
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapVisits'],0);
    %%
    figure
    VV=maps.(sen).dist.zonal.fromBirth.mean/1000;
    pcolor(lo,la,VV);shading flat
    caxis([ticks.dist(1) ticks.dist(2)])
    cb=decorate('dist',ticks,DD,sen,'Distance from Birth','km',0,1);
    doublemap(cb,win,aut,[1 1 1])
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapDFB'],0);
    %%
    figure
    VV=maps.(sen).dist.zonal.tillDeath.mean/1000;
    pcolor(lo,la,VV);shading flat
    caxis([ticks.dist(1) ticks.dist(2)])
    cb=decorate('dist',ticks,DD,sen,'Distance till Death','km',0,1);
    doublemap(cb,win,aut,[1 1 1])
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapDTD'],0);
    %%
    figure
    VV=log(maps.(sen).dist.traj.fromBirth.mean/1000);
    pcolor(lo,la,VV);shading flat
    caxis([ticks.disttot(1) ticks.disttot(2)])
    decorate('disttot',ticks,DD,sen,'Total distance travelled since birth','km',1,1);
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapDTFB'],0);
    %%
    figure
    VV=log(maps.(sen).dist.traj.tillDeath.mean/1000);
    pcolor(lo,la,VV);shading flat
    caxis([ticks.disttot(1) ticks.disttot(2)])
    decorate('disttot',ticks,DD,sen,'Total distance to be travelled till death','km',1,1);
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapDTTD'],0);
    %%
    figure
    VV=maps.(sen).vel.zonal.mean*100;
    pcolor(lo,la,VV);shading flat
    caxis([ticks.vel(1) ticks.vel(2)])
    cb=decorate('vel',ticks,DD,sen,'Zonal velocity','cm/s',0,1);
    doublemap(cb,aut,win,[.9 1 .9])
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapVel'],0);
    %%
    figure
    VV=maps.(sen).radius.mean.mean/1000;
    pcolor(lo,la,VV);shading flat
    caxis([ticks.radius(1) ticks.radius(2)])
    decorate('radius',ticks,DD,sen,'Radius','km',0,1);
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_MapRad'],0);
    
end

end
function trackPlots(DD,ticks,tracks)
%%
senses={'Cycs','AntiCycs'};
%figure
for sense=senses; sen=sense{1};
    
    drawColorLinem(tracks.(sen),'lat','isoper') ;
    title([sen '- deflections'])
    axis([-1200 300 -300 300])
    set(gca,'ytick',[-200 0 200])
    set(gca,'xtick',[-1000 0 200])
    colorbar
    xlabel('latitude repr. by color; IQ repr. by thickness')
    axis equal
    
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_defletcs']);
    
    
    
    %%
    field='age';
    drawColorLine(tracks.(sen),field,ticks.age(2),ticks.age(1),1,0) ;
    decorate(field,ticks,DD,sen,field,'d',1,1);
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_age']);
    %%
    %figure
    drawColorLine(tracks.(sen),'isoper',ticks.isoper(2),ticks.isoper(1),0,0) ;
    decorate('isoper',ticks,DD,sen,'IQ',' ',0,100)
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_IQ']);
    %%
    %figure
    drawColorLine(tracks.(sen),'radiusmean',ticks.radius(2)*1000,ticks.radius(1)*1000,0,0) ;
    decorate('radius',ticks,DD,sen,'Radius','km',0,1)
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_radius']);
    %%
    %figure
    
    drawColorLine(tracks.(sen),'peakampto_ellipse',ticks.amp(2)/100,ticks.amp(1)/100,0,0) ;
    decorate('amp',ticks,DD,sen,'Amp to ellipse','cm',1,1)
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_TrackPeakampto_ellipse']);
    %%
    %figure
    drawColorLine(tracks.(sen),'peakampto_contour',ticks.amp(2)/100,ticks.amp(1)/100,0,0) ;
    decorate('amp',ticks,DD,sen,'Amp to contour','cm',1,1)
    savefig(ticks.rez,ticks.width,ticks.height,[sen,'_TrackPeakampto_contour'])
end
end

function doublemap(cb,cm1,cm2,centercol)
%% get colorbardata
zlim=(get(cb,'ylim'));
ztick=(get(cb,'ytick'));
zticklabel=(get(cb,'yticklabel'));

%% resample to fit ticks
cm1r=resample(cm1,abs(zlim(1)),1);
cm2r=resample(cm2,zlim(2),1);
CM=[cm1r;flipud(cm2r)];
%% blend in the middle
midfilt=linspace(-1,zlim(2)/abs(zlim(1)),length(CM));
gp=repmat(gauspuls(midfilt,1,1,-1/5)',1,3);
centercolvec=repmat(centercol,size(CM,1),1);
CM=(1-gp).*CM + gp.*centercolvec;
%% correct for round errors
CM(CM<0)=0;
CM(CM>1)=1;
%% reset to old params
colormap(CM);
caxis(zlim);
set(cb,'ytick',ztick);
set(cb,'yticklabel',zticklabel);
end
function drawcoast
load coast;
hold on; plot(long,lat,'LineWidth',0.5);
end
function cb=decorate(field,ticks,DD,tit,tit2,unit,logornot,decim,coast)
if nargin<9
    coast=true;
end
axis(ticks.axis);
set(gca,'ytick',ticks.y);
set(gca,'xtick',ticks.x) ;
cb=colorbar;
if logornot
    zticks=linspace(log(ticks.(field)(1)),log(ticks.(field)(2)),ticks.(field)(3))';
    zticklabel=num2str(round(exp(zticks)));
else
    zticks=linspace(ticks.(field)(1),ticks.(field)(2),ticks.(field)(3))';
    zticklabel=num2str(round(zticks*decim)/decim);
end
caxis([zticks(1) zticks(end)])
set(cb,'ytick',zticks);
set(cb,'yticklabel',zticklabel);
title([tit,' - ',tit2,' [',unit,']'])
xlabel(['Eddies that died younger ',num2str(DD.thresh.life),' days are excluded'])
if coast
    drawcoast;
end
end




function [maxV,cmap]=drawColorLinem(files,fieldName,fieldName2)
cmap=jet;% Generate range of color indices that map to cmap

%% get extremata
minVglob=inf; minIQglob=inf; maxVglob=0; maxIQglob=0;
for file=files; ff=file{1};
    V=load(ff,fieldName,fieldName2);
    minV=min(V.(fieldName));
    maxV=max(V.(fieldName));
    minIQ=min(V.(fieldName2));
    maxIQ=max(V.(fieldName2));
    if minV<minVglob, minVglob=minV;end
    if minIQ<minIQglob, minIQglob=minIQ;end
    if maxV>maxVglob, maxVglob=maxV;end
    if maxIQ>maxIQglob, maxIQglob=maxIQ;end
end
iqiq=linspace(minIQglob,maxIQglob,10);
kk=linspace(minVglob,maxVglob,size(cmap,1));

for ee=1:numel(files)
    V=load(files{ee},fieldName,fieldName2,'lat','lon');
    VV=V.(fieldName);
    VViq=V.(fieldName2);
    if isempty(VV)
        continue
    end
    cm = spline(kk,cmap',VV);       % Find interpolated colorvalue
    cm(cm>1)=1;                        % Sometimes iterpolation gives values that are out of [0,1] range...
    cm(cm<0)=0;
    
    iq = spline(iqiq,linspace(0.01,1.5,10),VViq);       % Find interpolated thickness
    
    lo=V.lon;
    la=V.lat;
    %% deg2km
    yy=[0 cumsum(deg2km(diff(la)))];
    xx=[0 cumsum(deg2km(diff(lo)).*cosd((la(1:end-1)+la(2:end))/2))];
    for ii=1:length(xx)-1
        if  abs(xx(ii+1)-xx(ii))<1000 % avoid 0->360 jumps
            line([xx(ii) xx(ii+1)],[yy(ii) yy(ii+1)],'color',cm(:,ii),'LineWidth',iq(ii));
        end
    end
end
caxis([minV maxV])
end

function [maxV,cmap]=drawColorLine(files,fieldName,maxV,minV,logornot,zeroshift)
cmap=jet;% Generate range of color indices that map to cmap
if logornot
    maxV=log(maxV);
    minV(minV==0)=1;
    minV=log(minV);
end
kk=linspace(minV,maxV,size(cmap,1));
for ee=1:numel(files)
    V=load(files{ee},fieldName,'lat','lon');
    VV=V.(fieldName);
    if isempty(VV)
        continue
    end
    if logornot
        VV(VV==0)=1;
        cm = spline(kk,cmap',log(VV));  % Find interpolated colorvalue
    else
        cm = spline(kk,cmap',VV);       % Find interpolated colorvalue
    end
    cm(cm>1)=1;                        % Sometimes iterpolation gives values that are out of [0,1] range...
    cm(cm<0)=0;
    lo=V.lon;
    la=V.lat;
    if zeroshift
        lo=lo-lo(1);
        la=la-la(1);
    end
    
    for ii=1:length(la)-1
        if  abs(lo(ii+1)-lo(ii))<10 % avoid 0->360 jumps
            line([lo(ii) lo(ii+1)],[la(ii) la(ii+1)],'color',cm(:,ii),'LineWidth',0.8);
        end
    end
end
caxis([minV maxV])
end
