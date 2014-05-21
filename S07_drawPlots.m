%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 19-Apr-2014 17:39:11
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all figs are saved to ~/FIGS/  !
function S07_drawPlots
    %% init
    init_threads(8);
     spmd(8)
        [DD,threadData]=inits;
        %%	set ticks here!
        ticks.rez=100;
        ticks.width=297/25.4*ticks.rez;
        ticks.height=ticks.width * DD.dim.Y/DD.dim.X;        
%         ticks.height=ticks.width/sqrt(2); % Din a4
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
        ticks.minMax=cell2mat(extractfield( load([DD.path.analyzed.name, 'vecs.mat']), 'minMax'));
        %% main
        main(DD,threadData,ticks);
     end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,OUT]=inits
    DD=initialise();
    OUT.maps=load([DD.path.analyzed.name, 'maps.mat']);
    OUT.la=OUT.maps.Cycs.lat;
    OUT.lo=OUT.maps.Cycs.lon;
    %% collect tracks
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
    %% get vectors
    OUT.vecs=load([DD.path.analyzed.name, 'vecs.mat']);
end
%
%
function main(DD,IN,ticks)
    if labindex==1
        
        histstuff(IN.vecs,DD,ticks)
    end
    if labindex==2
        mapstuff(IN.maps,IN.vecs,DD,ticks,IN.lo,IN.la)
    end
    if labindex > 2
        trackPlots(DD,ticks,IN.tracks)
    end
end
%
%
function histstuff(vecs,DD,ticks)
    senses={'Cycs','AntiCycs'};
    for sense=senses;sen=sense{1};
        %%
        if isempty(vecs.(sen).lat), warning(['warning, no ' sen ' found!']);return;end 
        range.(sen).lat=round(nanmin(vecs.(sen).lat)):round(nanmax(vecs.(sen).lat));
        %%
        range.(sen).age=round(nanmin(vecs.(sen).age)):round(nanmax(vecs.(sen).age));
        %%
        range.(sen).cum=range.(sen).age;
        %%
        [hc.(sen).lat]=numPerLat(vecs.(sen).lat,DD,ticks,range.(sen).lat,sen);
        %%
        [hc.(sen).age,hc.(sen).cum]=ageCum(vecs.(sen).age,DD,ticks,range.(sen).age,sen);
    end
         
    %%
    ratioBar(hc,'lat',range,ticks);
    tit='ratio of cyclones to anticyclones as function of latitude';
    title(tit)
   savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,['RatLat'],0);
    %%
    ratioBar(hc,'age',range,ticks);
    tit='ratio of cyclones to anticyclones as function of age';
    xlab='age [d] - gray scale indicates total number of eddies available';
    title(tit);
    xlabel(xlab)   ;
   savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_LatNum'],0);
    %%
    ratioBar(hc,'cum',range,ticks);
    tit='ratio of cyclones to anticyclones as upper tail cum. of age';
    title(tit);
   savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,['RatCumAge'],0);
    
end
%
%
function mapstuff(maps,vecs,DD,ticks,lo,la)
    aut=autumn;
    win=winter;
    senses={'Cycs','AntiCycs'};
    for sense=senses;sen=sense{1};
        if isempty(vecs.(sen).lat), warning(['warning, no ' sen ' found!']);return;end 
        %%
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
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_deathsBirths'],0);
        
        %%
        figure
        maps.(sen).age.logmean=log(maps.(sen).age.mean);
        VV=maps.(sen).age.logmean;
        pcolor(lo,la,VV);shading flat
        caxis([ticks.age(1) ticks.age(2)])
        decorate('age',ticks,DD,sen,'age','d',1,1);
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_MapAge'],0)
        %%
        figure
        VV=maps.(sen).visits.single;
        VV(VV==0)=nan;
        pcolor(lo,la,log(VV));shading flat
        caxis([ticks.visits(1) ticks.visits(2)])
        decorate('visits',ticks,DD,sen,'Visits of unique eddy',' ',0,1);
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_MapVisits'],0);
        %%
        figure
        VV=maps.(sen).dist.zonal.fromBirth.mean/1000;
        pcolor(lo,la,VV);shading flat
        caxis([ticks.dist(1) ticks.dist(2)])
        cb=decorate('dist',ticks,DD,sen,'Distance from Birth','km',0,1);
         doublemap(cb,win,aut,[.9 1 .9])
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_MapDFB'],0);
        %%
        figure
        VV=maps.(sen).dist.zonal.tillDeath.mean/1000;
        pcolor(lo,la,VV);shading flat
        caxis([ticks.dist(1) ticks.dist(2)])
        cb=decorate('dist',ticks,DD,sen,'Distance till Death','km',0,1);
        doublemap(cb,win,aut,[.9 1 .9])
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_MapDTD'],0);
        %%
        figure
        VV=log(maps.(sen).dist.traj.fromBirth.mean/1000);
        pcolor(lo,la,VV);shading flat
        caxis([ticks.disttot(1) ticks.disttot(2)])
        decorate('disttot',ticks,DD,sen,'Total distance travelled since birth','km',0,1);
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_MapDTFB'],0);
        %%
        figure
        VV=log(maps.(sen).dist.traj.tillDeath.mean/1000);
        pcolor(lo,la,VV);shading flat
        caxis([ticks.disttot(1) ticks.disttot(2)])
        decorate('disttot',ticks,DD,sen,'Total distance to be travelled till death','km',0,1);
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_MapDTTD'],0);
        %%
        figure
        VV=maps.(sen).vel.zonal.mean*100;
        pcolor(lo,la,VV);shading flat
        caxis([ticks.vel(1) ticks.vel(2)])
        cb=decorate('vel',ticks,DD,sen,'Zonal velocity','cm/s',0,1);
        doublemap(cb,aut,win,[.9 1 .9])
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_MapVel'],0);
        %%
        figure
        VV=maps.(sen).radius.mean.mean/1000;
        pcolor(lo,la,VV);shading flat
        caxis([ticks.radius(1) ticks.radius(2)])
        decorate('radius',ticks,DD,sen,'Radius','km',0,1);
       savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_MapRad'],0);
        
    end
    
end
%
%
function trackPlots(DD,ticks,tracks)
    %%
    senses=fieldnames(tracks);
    for sense=senses; sen=sense{1};
        if labindex==3
            %%
            drawColorLinem(ticks,tracks.(sen),'lat','isoper') ;
            title([sen '- deflections'])
            axis([-1200 300 -300 300])
            set(gca,'ytick',[-200 0 200])
            set(gca,'xtick',[-1000 0 200])
            colorbar
            xlabel('latitude repr. by color; IQ repr. by thickness')
            axis equal
           savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_defletcs']);
        elseif labindex==4
            %%
            field='age';
            drawColorLine(ticks,tracks.(sen),field,ticks.age(2),ticks.age(1),1,0) ;
            decorate(field,ticks,DD,sen,field,'d',1,1);
           savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_age']);
        elseif labindex==5
            %%
            drawColorLine(ticks,tracks.(sen),'isoper',ticks.isoper(2),ticks.isoper(1),0,0) ;
            decorate('isoper',ticks,DD,sen,'IQ',' ',0,100)
           savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_IQ']);
        elseif labindex==6
            %%
            drawColorLine(ticks,tracks.(sen),'radiusmean',ticks.radius(2)*1000,ticks.radius(1)*1000,0,0) ;
            decorate('radius',ticks,DD,sen,'Radius','km',0,1)
           savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_radius']);
        elseif labindex==7
            %%
            drawColorLine(ticks,tracks.(sen),'peakampto_ellipse',ticks.amp(2)/100,ticks.amp(1)/100,0,0) ;
            decorate('amp',ticks,DD,sen,'Amp to ellipse','cm',1,1)
           savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_TrackPeakampto_ellipse']);
        elseif labindex==1
            %%
            drawColorLine(ticks,tracks.(sen),'peakampto_contour',ticks.amp(2)/100,ticks.amp(1)/100,0,0) ;
            decorate('amp',ticks,DD,sen,'Amp to contour','cm',1,1)
           savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_TrackPeakampto_contour'])
        end
    end
end
%
%
function logyBar(totnum,x,y)
    y=log10(y);
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
    cm=flipud(bone);
    colormap(cm)
    cb=colorbar;
    zticks=linspace(2,max(totnum),5)';
    zticklabel=num2str(round((zticks)));
    caxis([zticks(1) zticks(end)])
    set(cb,'ytick',zticks);
    set(cb,'yticklabel',zticklabel);
end
%
%
function [lat]=numPerLat(latin,DD,ticks,range,sen)
    figure
    lat=histc(latin,range);
    semilogy(range,lat);
    tit=['number of ',sen,' per 1° lat'];
    title([tit]);
   savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_latNum'],0);
end
%
%
function [age,cum]=ageCum(agein,DD,ticks,range,sen)
    figure
    age=histc(agein,range);
    semilogy(range,age)
    tit=['number of ',sen,' per age'];
    unit='d';
    title([tit,' [',unit,']'])
    cum=fliplr(cumsum(fliplr(age)));
    semilogy(range,cum)
    tit=['upper tail cumulative of ',sen,' per age'];
    title([tit,' [',unit,']'])
   savefig(DD.path.plots,ticks.rez,ticks.width,ticks.height,[sen,'_ageUpCum'],0);
end
%
%
function [totnum]=ratioBar(hc,field,range,ticks)
    %% ratio cyc/acyc per lat
    figure
    %% find max length for ratio vecotr
    len=nanmin([length(hc.Cycs.(field)),length(hc.AntiCycs.(field))]);
    hc.rat.(field)=hc.Cycs.(field)(1:len)./hc.AntiCycs.(field)(1:len);
    hc.rat.(field)(isinf(hc.rat.(field)))=nan;
    hc.rat.(field)((hc.rat.(field)==0))=nan;
    %% total length
    totnum=hc.Cycs.(field)(1:len) + hc.AntiCycs.(field)(1:len);
    %%%%
    logyBar(totnum,range.Cycs.(field)(1:len),hc.rat.(field))
    %%%%
    axis([range.Cycs.(field)(1) range.Cycs.(field)(end) -1 1]);
    ticks.y=linspace(-1,1,3);
    ticks.lab.y=num2cell(10.^(ticks.y));
    set(gca,'ytick',ticks.y)
    set(gca,'yticklabel',ticks.lab.y)
    set(gca,'xtick',linspFromTriple(ticks.lat))
end
%
%
function doublemap(cb,cm1,cm2,centercol)
    %% get colorbardata
    zlim=(get(cb,'ylim'));
    ztick=(get(cb,'ytick'));
    zticklabel=(get(cb,'yticklabel'));
   
    %% resample to fit ticks
    cm1r=resample(cm1,abs(zlim(1)),round(abs(zlim(1)) * abs(zlim(2)/zlim(1))));
    cm2r=resample(cm2,zlim(2),zlim(2));
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
%
%
function drawcoast
    load coast;
    hold on; plot(long,lat,'LineWidth',0.5);
end
%
%
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
%
%
function [maxV,cmap]=drawColorLinem(ticks,files,fieldName,fieldName2)
    cmap=jet;% Generate range of color indices that map to cmap
    %% get extremata
    maxIQ=ticks.isoper(2);
    minIQ=ticks.isoper(1);
    minV=ticks.lat(1);
    maxV=ticks.lat(2);
    iqiq=linspace(minIQ,maxIQ,10);
    kk=linspace(minV,maxV,size(cmap,1));
    for ee=1:numel(files)
        V=load(files{ee},fieldName,fieldName2,'lat','lon');
        VV=V.(fieldName);
        VViq=V.(fieldName2);
        if isempty(VV)
            continue
        end
        cm = spline(kk,cmap',VV);       % Find interpolated colorvalue
        cm(cm>1)=1;                     % Sometimes iterpolation gives values that are out of [0,1] range...
        cm(cm<0)=0;
        %% Find interpolated thickness
        iq = spline(iqiq,linspace(0.001,2,10),VViq);
        iq(iq>1)=1;                       
        iq(iq<0)=0.0000001;
        %% deg2km
        yy=[0 cumsum(deg2km(diff(V.lat)))];
        xx=[0 cumsum(deg2km(diff(V.lon)).*cosd((V.lat(1:end-1)+V.lat(2:end))/2))];
        for ii=1:length(xx)-1
            if  abs(xx(ii+1)-xx(ii))<1000 % avoid 0->360 jumps
                line([xx(ii) xx(ii+1)],[yy(ii) yy(ii+1)],'color',cm(:,ii),'LineWidth',iq(ii));
            end
        end
    end
    caxis([minV maxV])
end
%
%
function [maxV,cmap]=drawColorLine(ticks,files,fieldName,maxV,minV,logornot,zeroshift)
    cmap=jet;% Generate range of color indices that map to cmap
    if logornot
        maxV=log(maxV);
        minV(minV==0)=1;
        minV=log(minV);
    end
    kk=linspace(minV,maxV,size(cmap,1));
    %%
    maxIQ=ticks.isoper(2);
    minIQ=ticks.isoper(1);
    iqiq=linspace(minIQ,maxIQ,10);
    %%
    for ee=1:numel(files)
        V=load(files{ee},fieldName,'isoper','lat','lon');
        VV=V.(fieldName);
        VViq=V.isoper;
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
        %% Find interpolated thickness
        iq = spline(iqiq,linspace(0.001,2,10),VViq);
        iq(iq>1)=1;                       
        iq(iq<0)=0.0000001;
        %%
        for ii=1:length(la)-1
            if  abs(lo(ii+1)-lo(ii))<10 % avoid 0->360 jumps
                line([lo(ii) lo(ii+1)],[la(ii) la(ii+1)],'color',cm(:,ii),'LineWidth',iq(ii));
            end
        end
    end
    caxis([minV maxV])
end


