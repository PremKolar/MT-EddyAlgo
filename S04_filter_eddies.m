%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% walks through all the contours and decides whether they qualify
function S04_filter_eddies
    %% init
    DD=initialise('conts',mfilename);
    DD.threads.num=init_threads(DD.threads.num);
    rossby=getRossbyPhaseSpeedAndRadius(DD);
    %% spmd
    main(DD,rossby);
    %% update infofile
    conclude(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,rossby)
    if DD.debugmode
        spmd_body(DD,rossby)
    else
        spmd(DD.threads.num)
            spmd_body(DD,rossby)
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD,rossby)
    [JJ]=SetThreadVar(DD);
    Td=disp_progress('init','filtering contours');
    for jj=1:numel(JJ)
        Td=disp_progress('disp',Td,numel(JJ));
        [EE,skip]=work_day(DD,JJ(jj),rossby);
        if skip,disp(['skipping ' num2str(jj)]);return;end
        %% save
        save_eddies(EE);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EE,skip]=work_day(DD,JJ,rossby)
    %% check for exisiting data
    skip=false;
    EE.filename.cont=JJ.files;
    EE.filename.cut=[DD.path.cuts.name, DD.pattern.prefix.cuts, JJ.protos];
    EE.filename.self=[DD.path.eddies.name, DD.pattern.prefix.eddies ,JJ.protos];
    if exist(EE.filename.self,'file') && ~DD.overwrite, skip=true; return; end
    %% get ssh data
    try
        cut=load(EE.filename.cut);
        %% get contours
        cont=load(EE.filename.cont);
    catch failed
        fprintf('cannot read %s! \n',EE.filename.cont)
        system(['rm ' EE.filename.cont])
        disp(failed.message);
        save(sprintf('S04fail-%s.mat',datestr(now,'mmddHHMM')));
        skip = true; return;
    end
    %% put all eddies into a struct: ee(number of eddies).characteristica
    ee=eddies2struct(cont.all,DD.thresh.corners);
    %% remember date
    [ee(:).daynum]=deal(JJ.daynums);
    %% avoid out of bounds integer coordinates close to boundaries
    [ee_clean,cut]=CleanEDDies(ee,cut,DD.contour.step);
    %% find them
    EE=find_eddies(EE,ee_clean,rossby,cut,DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EE=find_eddies(EE,ee,rossby,cut,DD)
    %% anti cyclones
    senN=[-1 1];
    for ii=1:2
        sen=DD.FieldKeys.senses{ii};
        [EE.(sen),EE.pass.(sen)]=walkThroughContsVertically(ee,rossby,cut,DD,senN(ii));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eddies, pass]=walkThroughContsVertically(ee,rossby,cut,DD,sense)
    pp=0;  pass=initPass(numel(ee))    ;
    %% init
    [eddyType,Zloop]=determineSense(DD.FieldKeys.senses,sense,numel(ee));
    %% loop
    Tv=disp_progress('init','running through contours vertically');
    for kk=Zloop % dir dep. on sense
        Tv= disp_progress('disp',Tv,numel(Zloop),5);
        [pass(kk),ee_out]=run_eddy_checks(pass(kk),ee(kk),rossby,cut,DD,sense);
        if all(struct2array(pass(kk))), pp=pp+1;
            %% append healthy found eddy TODO
            eddies(pp)=ee_out; %#ok<*AGROW>
            %% flag respective overlap too
            if strcmp(cut.window.type,'globe')
                [yi,xi]=find(ee_out.mask);
                doubleFlag.east=xi>cut.window.size.X;
                doubleFlag.west=xi<=cut.window.sizePlus.X - cut.window.size.X  +1;
                xi = [ xi xi(doubleFlag.east)-cut.window.size.X xi(doubleFlag.west)+cut.window.size.X ];
                yi = [ yi yi(doubleFlag.east)                   yi(doubleFlag.west)                   ];
                ee_out.mask(drop_2d_to_1d(yi,xi,cut.window.size.Y))=true;
            end
            %% nan out ssh where eddy was found
            cut.grids.ssh(ee_out.mask)=nan;
        end
    end
    %% catch
    if ~any(struct2array(pass(:)))
        error('no %s made it through the filter...',eddyType)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eddyType,Zloop]=determineSense(senseKeys,sense,NumEds)
    switch sense
        case -1
            eddyType=senseKeys(1); % anti cycs
            Zloop=1:NumEds;
        case 1
            eddyType=senseKeys(2); %  cycs
            Zloop=NumEds:-1:1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass,ee]=run_eddy_checks(pass,ee,rossby,cut,DD,direction)
    %% pre-nan-check
    pass.rim=CR_RimNan(ee.coordinates.int, cut.dim.Y, cut.grids.ssh);
    if ~pass.rim, return, end;
    %% closed ring check
    [pass.CR_ClosedRing]=CR_ClosedRing(ee);
    if ~pass.CR_ClosedRing, return, end;
    %% pre filter 'thin 1dimensional' eddies
    pass.CR_2dEDDy=CR_2dEDDy(ee.coordinates.int);
    if ~pass.CR_2dEDDy, return, end;
    %% get coordinates for zoom cut
    [zoom,pass.winlim]=get_window_limits(ee.coordinates,4,DD.map.window);
    if ~pass.winlim, return, end;
    %% cut out rectangle encompassing eddy range only for further calcs
    zoom.fields=EDDyCut_init(cut.grids,zoom);
    %% generate logical masks defining eddy interiour and outline
    zoom.mask=EDDyCut_mask(zoom);
    %% check for nans matlab.matwithin eddy
    [pass.CR_Nan]=CR_Nan(zoom);
    if ~pass.CR_Nan, return, end;
    if ~any(zoom.mask.inside),pass.CR_Nan=false; return; end
    %% check for correct sense
    [pass.CR_sense,ee.sense]=CR_sense(zoom,direction,ee.level);
    if ~pass.CR_sense, return, end;
    %% calculate area with respect to contour
    RoL=getLocalRossyRadius(rossby.Lr,ee.coordinates.int);
    [ee.area,pass.Area]=Area(zoom,RoL,DD.thresh.maxRadiusOverRossbyL);
    if ~pass.Area && DD.switchs.maxRadiusOverRossbyL, return, end;
    %% calc contour circumference in [SI]
    [ee.circum.si,ee.fourierCont]=EDDyCircumference(zoom);
    %% filter eddies not circle-like enough
    [pass.CR_Shape,ee.isoper, ee.chelt]=CR_Shape(zoom,ee,DD.thresh.shape,DD.switchs);
    if ~pass.CR_Shape, return, end;
    %% get peak position and amplitude w.r.t contour
    [pass.CR_AmpPeak,ee.peak,zoom.ssh_BasePos]=CR_AmpPeak(ee,zoom,DD.thresh.amp);
    if ~pass.CR_AmpPeak, return, end;
    %% get profiles
    [ee.profiles,~]=EDDyProfiles(ee,zoom);
    %% get radius according to max UV ie min vort
    [ee.radius,pass.CR_radius]=EDDyRadiusFromUV(ee.peak.z, ee.profiles,DD.thresh.radius);
    
    %     if a
    %         [ee.profiles,~]=EDyrp(ee,zoom);
    %     end
    %
    if ~pass.CR_radius, return, end;
    %% get ideal ellipse contour
    zoom.mask.ellipse=EDDyEllipse(ee,zoom.mask);
    %% get effective amplitude relative to ellipse;
    ee.peak.amp.to_ellipse=EDDyAmp2Ellipse(ee,zoom);
    %% append mask to ee in cut coordinates
    [ee.mask]=sparse(EDDyPackMask(zoom.mask.filled,zoom.limits,cut.dim));
    %% get center of 'volume'
    [ee.volume]=CenterOfVolume(zoom,ee.area.total,cut.dim.Y);
    %% get area centroid (chelton style)
    [ee.centroid]=AreaCentroid(zoom,cut.dim.Y);
    %% get coordinates
    [ee.geo]=geocoor(zoom,ee.volume);
    %% append 'age'
    ee.age=0;
    %% get trackref
    ee.trackref=getTrackRef(ee,DD.parameters.trackingRef);
    %% append projected location
    if (DD.switchs.distlimit && DD.switchs.RossbyStuff)
        [ee.projLocsMask]=ProjectedLocations(rossby.c,cut,DD,ee.trackref);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RoL=getLocalRossyRadius(rossbyL,coor)
    [Y,X]=size(rossbyL);
    x=round(coor.x);
    y=round(coor.y);
    x(x<1)=1;
    x(x>X)=X;
    y(y<1)=1;
    y(y>Y)=Y;
    RoL=nanmedian(rossbyL(drop_2d_to_1d(y,x,Y)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass,sense]=CR_sense(zoom,direc,level)
    pass=false;
    sense=struct;
    %% water column up: seeking anti cyclones; down: cyclones
    if direc==-1
        if all(zoom.fields.ssh(zoom.mask.inside) >= level )
            pass=true;
            sense.str='AntiCyclonic';
            sense.num=-1;
        end
    elseif direc==1
        if all(zoom.fields.ssh(zoom.mask.inside) <= level )
            pass=true;
            sense.str='Cyclonic';
            sense.num=1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pass=CR_RimNan(coor, Y, ssh)
    pass=true;
    if any(isnan(ssh(drop_2d_to_1d(coor.y, coor.x, Y)))), pass=false; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass,peak,base]=CR_AmpPeak(ee,z,thresh)
    pass=false;
    peak.mean_ssh=mean(z.fields.ssh(z.mask.filled));
    %% make current level zero level and zero out everything else
    base=poslin(-ee.sense.num*(z.fields.ssh-ee.level));
    base(~z.mask.filled)=0;
    %% amplitude
    [peak.amp.to_contour,peak.lin]=max(base(:));
    [peak.z.y,peak.z.x]=raise_1d_to_2d(z.size.Y, peak.lin);
    peak.amp.to_mean = z.fields.ssh(peak.lin)-peak.mean_ssh;
    %% coordinates in full map
    peak.y=peak.z.y+z.limits.y(1) -1;
    peak.x=peak.z.x+z.limits.x(1) -1;
    if peak.amp.to_contour>=thresh,	pass=true; 	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass,IQ,chelt]=CR_Shape(z,ee,thresh,switches)
    [passes.iq,IQ]=IsopQuo(ee,thresh.iq);
    [passes.chelt,chelt]=chelton_shape(z,ee);
    if switches.IQ && ~switches.chelt
        pass=passes.iq;
    elseif switches.chelt && ~switches.IQ
        pass=passes.chelt;
    elseif switches.chelt && switches.IQ
        pass=passes.chelt && passes.iq;
    else
        error('choose at least one shape method (IQ or chelton method in input_vars switches section)') %#ok<*ERTAG>
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass,chelt]=chelton_shape(z,ee)
    %% get max dist(all2all)
    f=ee.fourierCont;
    x=f.x(f.ii);
    y=f.y(f.ii);
    xiy=x+1i*y;
    [A,B]=meshgrid(xiy,xiy);
    maxDist=max(max(abs(A-B)))*1e3;
    %%
    medlat=abs(nanmean(reshape(z.mask.rim_only.*z.fields.lat,1,[]))) ;
    %%
    if medlat> 25
        chelt  = 1 - maxDist/4e5;
    else
        chelt  =  1 - maxDist/(8e5*(25 - medlat)/25 + 4e5);
    end
    if chelt >= 0, pass=true; else pass=false; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass,iq]=IsopQuo(ee,thresh)
    %% isoperimetric quotient
    getIQ  = @(area,circum) 4*pi*area/circum^2;
    iq = getIQ(ee.area.intrp,ee.circum.si);
    %%
    if iq >= thresh, pass=true; else pass=false; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass]=CR_2dEDDy(coor)
    if (max(coor.x)-min(coor.x)<2) || (max(coor.y)-min(coor.y)<2)
        pass=false;
    else
        pass=true;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass]=CR_Nan(z)
    ssh=z.fields.ssh(z.mask.filled);
    if ~any(isnan(ssh(:))), pass=true; else pass=false; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass]=CR_ClosedRing(ee)
    x=ee.coordinates.int.x;
    y=ee.coordinates.int.y;
    if abs(x(1)-x(end))>1 || abs(y(1)-y(end))>1;
        pass=false;
    else
        pass=true;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RS=getRossbyPhaseSpeedAndRadius(DD)
    if DD.switchs.RossbyStuff
        RS.Lr=getfield(load([DD.path.Rossby.name 'RossbyRadius.mat']),'data');
        RS.c=getfield(load([DD.path.Rossby.name 'RossbyPhaseSpeed.mat']),'data');
    else
        warning('No Rossby Radius available. Ignoring upper constraint on eddy scale!') %#ok<*WNTAG>
        RS.c=[];
        RS.Lr=[];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [centroid]=AreaCentroid(zoom,Y)
    %% factor each grid cell equally (compare to CenterOfVolume())
    ssh=double(logical(zoom.ssh_BasePos));
    %% get centroid:   COVs = \frac{1}{A} \sum_{i=1}^n 1 \vec{x}_i,
    [XI,YI]=meshgrid(1:size(ssh,2), 1:size(ssh,1));
    y=sum(nansum(ssh.*YI));
    x=sum(nansum(ssh.*XI));
    yz=(y/nansum(ssh(:)));
    xz=(x/nansum(ssh(:)));
    y=yz + double(zoom.limits.y(1))-1;
    x=xz + double(zoom.limits.x(1))-1;
    centroid.xz=xz;
    centroid.yz=yz;
    centroid.x=x;
    centroid.y=y;
    centroid.lin=drop_2d_to_1d(y,x,Y);
    centroid.linz=drop_2d_to_1d(yz,xz,size(ssh,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask,trackref]=ProjectedLocations(rossbyU,cut,DD,trackref)
    %% get rossby wave phase speed
    
    rU=rossbyU(trackref.lin);
    rU(abs(rU)>1)=sign(rU)*1;  % put upper limit on rossby wave speed
    %% get projected distance (1.75 * dt*rU  as in chelton 2011)
    dist.east=DD.parameters.minProjecDist;
    dist.southnorth=dist.east;
    dist.ro=abs(DD.time.delta_t*86400* DD.parameters.rossbySpeedFactor * rU);
    dist.west=max([dist.ro, dist.east]); % not less than dist.east, see chelton 11
    %% get major/minor semi-axes [m]
    ax.maj=sum([dist.east, dist.west])/2;
    ax.min=DD.parameters.minProjecDist;
    %% get dx/dy at that eddy pos
    dx=cut.grids.DX(trackref.lin);
    dy=cut.grids.DY(trackref.lin);
    %% get major/minor semi-axes [increments]
    ax.majinc=ceil(ax.maj/dx);
    ax.mininc=ceil(ax.min/dy);
    %% translate dist to increments
    dist.eastInc=(ceil(dist.east/dx));
    dist.westInc=(ceil(dist.west/dx));
    dist.southnorthInc=(ceil(dist.southnorth/dy));
    %% get positions of params
    xi.f2=(trackref.x);
    yi.f2=(trackref.y);
    xi.center=round(xi.f2-(ax.majinc-dist.eastInc));
    yi.center=round(yi.f2);
    %% build x vector (major axis >= minor axis always!)
    fullcirc=linspace(0,2*pi,4*numel(-dist.westInc:dist.eastInc));
    ellip.x=round(ax.majinc * cos(fullcirc)) + xi.center;
    ellip.y=round(ax.mininc * sin(fullcirc)) + yi.center;
    ellip.lin=unique(drop_2d_to_1d(ellip.y,ellip.x,cut.dim.Y));
    %% take care of out of bounds values (only applicable to zonally non continous case. this shouldnt happen in global case)
    ellip.x(ellip.x<1)=1;
    ellip.x(ellip.x>cut.dim.X)=cut.dim.X;
    ellip.y(ellip.y<1)=1;
    ellip.y(ellip.y>cut.dim.Y)=cut.dim.Y;
    xi.center(xi.center<1)=1;
    xi.center(xi.center>cut.dim.X)=cut.dim.X;
    yi.center(yi.center<1)=1;
    yi.center(yi.center>cut.dim.Y)=cut.dim.Y;
    %% flag respective overlap too
    if strcmp(cut.window.type,'globe')
        doubleFlag.east=ellip.x>cut.window.size.X;
        doubleFlag.west=ellip.x<=cut.window.sizePlus.X - cut.window.size.X  +1;
        ellip.x =[ ellip.x ellip.x(doubleFlag.east)-cut.window.size.X ellip.x(doubleFlag.west)+cut.window.size.X ];
        ellip.y =[ ellip.y ellip.y(doubleFlag.east)                   ellip.y(doubleFlag.west)                   ];
    end
    %% build boundary mask
    mask.logical=false(struct2array(cut.dim));
    mask.logical(drop_2d_to_1d(ellip.y,ellip.x,cut.dim.Y))=true;
    mask.logical=sparse(imfill(mask.logical,double([yi.center xi.center]),4));
    mask.lin=find(mask.logical);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TR=getTrackRef(ee,tr)
    switch tr
        case 'centroid'
            TR=ee.centroid;
        case 'CenterOfVolume'
            TR=ee.volume.center;
        case 'peak'
            TR=ee.peak;
    end
    
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_eddies(EE)
    save(EE.filename.self,'-struct','EE')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [area,pass]=Area(z,rossbyL,scaleThresh)
    area=struct;
    area.pixels=(z.fields.DX.*z.fields.DY).*(z.mask.inside + z.mask.rim_only/2);  % include 'half of rim'
    area.total=sum(area.pixels(:));
    area.meanPerSquare=mean(z.fields.DX(z.mask.filled).*z.fields.DY(z.mask.filled));
    area.intrp=area.meanPerSquare*polyarea(z.coor.exact.x,z.coor.exact.y);
    area.RadiusOverRossbyL=sqrt(area.intrp/pi)/rossbyL;
    if area.RadiusOverRossbyL > scaleThresh
        pass=false;
    else
        pass=true;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function	[mask_out]=EDDyPackMask(mask_in,limits,dims)
    mask_out=false(dims.Y,dims.X);
    mask_out(limits.y(1):limits.y(2),limits.x(1):limits.x(2))=mask_in;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [amp] = EDDyAmp2Ellipse(ee,zoom)
    %% mean amplitude with respect to ellipse contour
    halfstep=@(data,y,x,dy,dx) mean([data(y,x) data(y+dy*2,x+dx*2)]);
    xa=double(ee.radius.coor.Xwest);
    xb=double(ee.radius.coor.Xeast);
    ya=double(ee.radius.coor.Ysouth);
    yb=double(ee.radius.coor.Ynorth);
    cx=double(ee.peak.z.x);
    cy=double(ee.peak.z.y);
    clear ssh
    ssh.west=halfstep(zoom.fields.ssh,cy,xa,0,.5);
    ssh.east=halfstep(zoom.fields.ssh,cy,xb,0,-.5);
    ssh.south=halfstep(zoom.fields.ssh,ya,cx,.5,0);
    ssh.north=halfstep(zoom.fields.ssh,yb,cx,-.5,0);
    ssh.mean=mean(struct2array(ssh));
    amp=abs(zoom.fields.ssh(ee.peak.z.y,ee.peak.z.x)- ssh.mean);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ellipse]=EDDyEllipse(ee,mask)
    %% get center, minor and major axis for ellipse
    xa=ee.radius.coor.Xwest;
    xb=ee.radius.coor.Xeast;
    ya=ee.radius.coor.Ysouth;
    yb=ee.radius.coor.Ynorth;
    xm=(mean([xa,xb]));
    ym=(mean([ya,yb]));
    axisX=(double(xb-xa))/2;
    axisY=(double(yb-ya))/2;
    %% init ellipse mask
    ellipse=false(mask.size.Y,mask.size.X);
    %% get ellipse coordinates
    linsdeg=(linspace(0,2*pi,2*sum(struct2array(mask.size))));
    ellipseX=round(axisX*cos(linsdeg)+xm);
    ellipseY=round(axisY*sin(linsdeg)+ym);
    ellipseX(ellipseX>mask.size.X)=mask.size.X;
    ellipseY(ellipseY>mask.size.Y)=mask.size.Y;
    ellipseX(ellipseX<1)=1;
    ellipseY(ellipseY<1)=1;
    xlin=unique(drop_2d_to_1d(ellipseY,ellipseX,mask.size.Y));
    %% draw into mask
    ellipse(xlin)=true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outIdx]=avoidLand(ssh,p)
    land  = reshape(isnan([nan reshape(ssh,1,[]) nan]),1,[]);
    ii  = [1 1:numel(ssh) numel(ssh)];
    a = ii(find( ii <= p & land,1,'last')  +1) ;
    b = ii(find( ii >= p & land,1,'first') -1) ;
    outIdx=a:b;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s,f]=EDDyProfiles(ee,z)
    %% detect meridional and zonal profiles shifted to baselevel of current level
    offset_term=ee.peak.amp.to_contour*ee.sense.num-ee.level;
    %%	zonal cut
    ssh =  -ee.sense.num*(z.fields.ssh(ee.peak.z.y,:) + offset_term);
    [wtr.x] = avoidLand(ssh,ee.peak.z.x);
    prof.x.ssh = (ssh(wtr.x));
    prof.x.ddis=z.fields.DX(ee.peak.z.y,wtr.x) ;
    prof.x.dist=z.fields.km_x(ee.peak.z.y,wtr.x)*1e3 ;
    %% meridional cut
    ssh =-ee.sense.num * (z.fields.ssh(:,ee.peak.z.x) + offset_term);
    wtr.y = avoidLand(ssh,ee.peak.z.y);
    prof.y.ssh = (ssh(wtr.y));
    prof.y.ddis=z.fields.DY(wtr.y,ee.peak.z.x) ;
    prof.y.dist=z.fields.km_y(wtr.y,ee.peak.z.x)*1e3 ;
    %
    %%	cranck up res
    for xyc={'x','y'};xy=xyc{1};
        s.(xy).dist = linspace(prof.(xy).dist(1), prof.(xy).dist(end),1e3)';
        s.(xy).idx  = linspace(wtr.(xy)(1), wtr.(xy)(end),1e3)';
        f.(xy)		= spline(prof.(xy).dist',prof.(xy).ssh');
        s.(xy).ssh  = (ppval(f.(xy), s.(xy).dist));
    end
    %% intrp versions
    fs=@(diflevel,arg,val) diffCentered(diflevel,arg,val)';
    f.FFour.x.ssh = fit(s.x.dist, (s.x.ssh),'fourier4');
    f.FFour.y.ssh = fit(s.y.dist, (s.y.ssh), 'fourier4');
    f.FOne.x.ssh = fit(s.x.dist, (s.x.ssh),'fourier1'); % can be used for better detrending
    f.FOne.y.ssh = fit(s.y.dist, (s.y.ssh), 'fourier1');
    
    %%
    for xyc={'x','y'};		xy=xyc{1};
        for foc={'FFour'};		fo=foc{1};
            s.(fo).(xy).sshf = feval(f.(fo).(xy).ssh, s.(xy).dist);
            s.(fo).(xy).UV   = fs(1,s.(xy).dist,s.(fo).(xy).sshf);
            s.(fo).(xy).UVd  = fs(2,s.(xy).dist,s.(fo).(xy).sshf);
        end
    end
    
    %%
    %             [~,pex]=min(abs(s.x.idx - double(ee.peak.z.x)));
    %             [~,pey]=min(abs(s.y.idx - double(ee.peak.z.y)));%
    %             figure('visible','off')
    %
    %             subplot(211)
    %             hold on
    %                 plot(s.x.dist,s.x.ssh,s.x.dist,s.FEight.x.sshf,s.x.dist,s.FFour.x.sshf)
    %             plot(s.x.dist,s.x.ssh,s.x.dist,s.FFour.x.sshf)
    %             axis tight
    %             ax=axis;
    %             plot(s.x.dist([pex pex]),ax([3 4]))
    %             subplot(212)
    %             hold on
    %             plot(s.y.dist,s.y.ssh,s.y.dist,s.FEight.y.sshf,s.y.dist,s.FFour.y.sshf)
    %             plot(s.y.dist,s.y.ssh,s.y.dist,s.FFour.y.sshf)
    %             axis tight
    %             ax=axis;
    %             plot(s.y.dist([pey pey]),ax([3 4]))
    %         saveas(gcf,[num2str(now),'.png'])
    %
    %     close all
    %%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [radius,pass]=EDDyRadiusFromUV(peak,prof,thresh)
    %%
    [x]=findSigma(peak,prof,'x');
    [y]=findSigma(peak,prof,'y');
    %%
    radius.zonal      = diff(x.dis([x.sigma-.5]))/2;
    radius.meridional = diff(y.dis([y.sigma-.5]))/2;
    radius.mean       = mean([radius.zonal;radius.meridional]);
    %% coors on ori grid
    halfidx=@(idx,sig) round(mean(idx([sig-.5 sig+.5])));
    radius.coor.Xwest  = halfidx(x.idx,x.sigma(1));
    radius.coor.Xeast  = halfidx(x.idx,x.sigma(2));
    radius.coor.Ysouth = halfidx(y.idx,y.sigma(1));
    radius.coor.Ynorth = halfidx(y.idx,y.sigma(2));
    %%
    
    %
    %     clf
    %     nrmc= @(x) (x-min(x))/max(x-min(x));
    %     figure(10000);clf;set(gcf,'visible','off');
    %     subplot(211)
    %     nrmdssh=nrmc(x.ssh);
    %     plot(x.idx,nrmdssh); hold on
    %     plot(x.idx([x.peakHigh x.peakHigh]),[0 1])
    %     plot(x.idx([x.peakLow x.peakLow]),[0 1],'r')
    %     subplot(212)
    %     nrmdssh=nrmc(y.ssh);
    %     plot(y.idx,nrmdssh); hold on
    %     plot(y.idx([y.peakHigh y.peakHigh]),[0 1])
    %     plot(y.idx([y.peakLow y.peakLow]),[0 1],'r')
    %     saveas(gcf,['AA',num2str(now),'.png'])
    %%
    
    %
    %     figure(1000)
    %     xy='y';
    %     nrmc= @(x) (x-min(x))/max(x-min(x));
    %     spl=@(x,abl) spline(1:abl,x,linspace(1,abl,100));
    %     pl  = @(x,ab) plot(nrmc(spl(x(ab(1):ab(2)),diff(ab)+1)));
    %     subplot(131)
    %     pl(prof.FFour.(xy).sshf,y.sigma);hold on;	grid minor;axis off tight
    %     subplot(132)
    %     pl(prof.FFour.(xy).UV,y.sigma	);hold on;	grid minor;axis off tight
    %     subplot(133)
    %     pl(prof.FFour.(xy).UVd,y.sigma	);hold on;	grid minor;axis off tight
    %     figure(5000)
    %     subplot(131)
    %     plot(prof.FFour.(xy).sshf);grid minor;axis  tight
    %     subplot(132)
    %     plot(prof.FFour.(xy).UV	);grid minor;axis  tight
    %     subplot(133)
    %     plot(prof.FFour.(xy).UVd	);grid minor;axis  tight
    %     %
    %%
    if radius.mean >= thresh, pass=true; else pass=false; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A]=findSigma(peak,prof,yx)
    % note that all profiles are 'high pressure' regardless of sense or
    % hemisphere here
    signJump.east=@(X) logical(diff(sign(X([1:end end]))));
    signJump.west=@(X) logical(diff(sign(X([1   1:end]))));
    %% x dir
    A.dUVdx = prof.FFour.(yx).UVd ;
    A.UV= prof.FFour.(yx).UV  ;
    A.ssh  = prof.FFour.(yx).sshf;
    A.dis  = prof.(yx).dist   ;
    A.idx  = prof.(yx).idx    ;
    A.intrpLen = length(A.idx);
    [~,A.peakLow] = min(abs(A.idx - double(peak.(yx)))); % index of peak in high res coors
    idxL=(1:length(A.ssh))';
    switch sign(A.UV(A.peakLow))
        case 1
            A.peakHigh = find(signJump.east(A.UV) &  idxL >= A.peakLow ,1,'first');
        case -1
            A.peakHigh = find(signJump.west(A.UV) &  idxL <= A.peakLow ,1,'last' );
    end
    
    if isempty(A.peakHigh)
        ergerg
    end
    
    % either left bndry or idx left of peak where dVdx crosses x-axis and slope of SSH is uphill
    F.a = A.idx < peak.(yx);
    F.b = signJump.west(A.dUVdx);
    F.c =  A.UV > 0;
    A.sigma(1) = max([ 1  find(F.a & F.b & F.c, 1, 'last') ]) +.5;
    % respectively for downhill side
    F.a = A.idx > peak.(yx);
    F.b = signJump.east(A.dUVdx); % right side of idx
    F.c =  A.UV < 0;
    A.sigma(2) = min([ numel(A.idx) find(F.a & F.b & F.c, 1, 'first')]) -.5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [geo]=geocoor(zoom,volume)
    xz=volume.center.xz;
    yz=volume.center.yz;
    geo.lat=interp2(zoom.fields.lat,xz,yz);
    geo.lon=interp2(zoom.fields.lon,xz,yz);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [volume]=CenterOfVolume(zoom,area,Y)
    %% get "volume" of eddy
    ssh=zoom.ssh_BasePos;
    volume.total=mean(ssh(:))*area;
    %% get center of volume  formula:   COVs = \frac{1}{M} \sum_{i=1}^n m_i \vec{x}_i,
    [XI,YI]=meshgrid(1:size(ssh,2), 1:size(ssh,1));
    y=sum(nansum(ssh.*YI));
    x=sum(nansum(ssh.*XI));
    yz=(y/nansum(ssh(:)));
    xz=(x/nansum(ssh(:)));
    y=yz + double(zoom.limits.y(1))-1;
    x=xz + double(zoom.limits.x(1))-1;
    volume.center.xz=xz;
    volume.center.yz=yz;
    volume.center.x=x;
    volume.center.y=y;
    volume.center.lin=drop_2d_to_1d(y,x,Y);
    volume.center.linz=drop_2d_to_1d(yz,xz,size(ssh,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [circum,f]=EDDyCircumference(z)
    %%
    ilin =z.coor.int.lin;
    x = z.fields.km_x(ilin);
    y = z.fields.km_y(ilin);
    f.ii=linspace(0,2*pi,numel(x))';
    f.II=linspace(0,2*pi,360)';
    options = fitoptions('Method','Smooth','SmoothingParam',0.99);
    f.x = fit(f.ii,x,'smoothingspline',options);
    f.y = fit(f.ii,y,'smoothingspline',options);
    circum=sum(hypot(diff(feval(f.x,f.II)),diff(feval(f.y,f.II)))) * 1000;
    
    % 	%%
    % 	clf
    % 	figure(1)
    % 	subplot(121)
    % 	plot(x,y,'r')
    % 	hold on
    % 	plot(f.x(f.II),f.y(f.II))
    % 	subplot(122)
    % 	hold on
    % 	plot(f.II(2:end),diff(f.x(f.II))','b',f.II(2:end),diff(f.y(f.II)),'r')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask=EDDyCut_mask(zoom)
    %% init
    dummymask=false(size(zoom.fields.ssh));
    [queryY,queryX]=find(~dummymask);
    queryLin = drop_2d_to_1d(queryY,queryX,size(dummymask,1));
    rimIntLin= drop_2d_to_1d(zoom.coor.int.y,zoom.coor.int.x,size(dummymask,1));
    %% inside
    querypoints=[queryX,queryY];
    node=struct2array(zoom.coor.exact);
    insideLin = queryLin(inpoly(querypoints,node)); % MAIN BOTTLENECK!!!!!
    mask.inside = dummymask;
    mask.inside(insideLin)=true;
    %% on rim
    mask.rim_only = dummymask;
    mask.rim_only(rimIntLin) = true;
    %% full
    mask.filled = mask.rim_only | mask.inside;
    %% dims
    [mask.size.Y, mask.size.X]=size(dummymask);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fields_out=EDDyCut_init(fields_in,z)
    ya=z.limits.y(1);
    yb=z.limits.y(2);
    xa=z.limits.x(1);
    xb=z.limits.x(2);
    %% cut all fields
    for ff=fieldnames(fields_in)'
        field=ff{1};
        fields_out.(field)=fields_in.(field)(ya:yb,xa:xb);
    end
    %%
    fields_out.km_x=deg2km(fields_out.lon);
    fields_out.km_y=deg2km(fields_out.lat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,passout]=get_window_limits(coor,enlargeFac,map)
    pass=true(3,1);
    z.coor=coor;
    %% output
    z.limits.x(1)=min(coor.int.x);
    z.limits.y(1)=min(coor.int.y);
    z.limits.x(2)=max(coor.int.x);
    z.limits.y(2)=max(coor.int.y);
    [z.limits,z.M]=enlarge_window(z.limits,enlargeFac,map.sizePlus) ;
    
    %%
    z.size.X=diff(z.limits.x)+1;
    z.size.Y=diff(z.limits.y)+1;
    z.coor.int.x=z.coor.int.x-z.limits.x(1) +1;
    z.coor.int.y=z.coor.int.y-z.limits.y(1) +1;
    z.coor.int.lin=drop_2d_to_1d(z.coor.int.y,z.coor.int.x,z.size.Y)	;
    z.coor.exact.x=z.coor.exact.x -double(z.limits.x(1))  +1;
    z.coor.exact.y=z.coor.exact.y -double(z.limits.y(1))  +1;
    %%
    if strcmp(map.type,'globe')
        %% in global case dismiss eddies touching zonal boundaries (another copy of these eddies exists that is not touching boundaries, due to the zonal appendage in S00b
        pass(1) =  z.limits.x(1)~=1;
        pass(2) =  z.limits.x(2)~=map.sizePlus.X;
        %% also dismiss eddies the zoom window of which is entirely inside the appended stripe ie beyond X_real(2). another legitimate copy of these exists in the western part of the actual map.)
        %         pass(3) =  ~all(coor.int.x > map.size.X);
    end
    passout=all(pass);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [inout,M]=enlarge_window(inout,factor,dim)
    
    half_width =round((diff(inout.x)+1)*(factor-1)/2);
    half_height=round((diff(inout.y)+1)*(factor-1)/2);
    inout.x(1)=max([1 inout.x(1)-half_width]);
    inout.x(2)=min([dim.X inout.x(2)+half_width]);
    inout.y(1)=max([1 inout.y(1)-half_height]);
    inout.y(2)=min([dim.Y inout.y(2)+half_height]);
    [M.x,M.y]=meshgrid(inout.x(1):inout.x(end),inout.y(1):inout.y(end));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EE]=eddies2struct(CC,thresh)
    EE=struct;
    ii=1;cc=0;
    while ii<size(CC,1);
        len=  CC(ii,2);% contourc saves the length of each contour before appending the next
        if len>=thresh.min && len<=thresh.max
            cc=cc+1;
            EE(cc).level=CC(ii,1);
            EE(cc).circum.length= len;
            EE(cc).coordinates.exact.x=CC(1+ii:ii+EE(cc).circum.length,1);
            EE(cc).coordinates.exact.y=CC(1+ii:ii+EE(cc).circum.length,2);
            EE(cc).coordinates.int.x=int32(EE(cc).coordinates.exact.x);
            EE(cc).coordinates.int.y=int32(EE(cc).coordinates.exact.y);
        end
        ii=ii+len+1; % jump to next eddy for next iteration
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ee,cut]=CleanEDDies(ee,cut,contstep) %#ok<INUSD>
    [cut.dim.Y,cut.dim.X]=size(cut.grids.ssh);
    for jj=1:numel(ee)
        x=ee(jj).coordinates.int.x;
        y=ee(jj).coordinates.int.y;
        %%
        x(x>cut.dim.X)=cut.dim.X;
        y(y>cut.dim.Y)=cut.dim.Y;
        x(x<1)=1;
        y(y<1)=1;
        %%
        ee(jj).coordinates.int.x=x;
        ee(jj).coordinates.int.y=y;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  pass=initPass(len)
    pass(len).rim=0;
    pass(len).CR_ClosedRing=0;
    pass(len).CR_2dEDDy=0;
    pass(len).winlim=0;
    pass(len).CR_Nan=0;
    pass(len).CR_sense=0;
    pass(len).Area=0;
    pass(len).CR_Shape=0;
    pass(len).CR_AmpPeak=0;
    pass(len).CR_radius=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
