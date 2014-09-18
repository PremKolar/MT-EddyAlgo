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
    dbstop if error
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
        
        toBeTried(DD,rossby,JJ,jj);
        
        %         try
        %             toBeTried(DD,rossby,JJ,jj);
        %         catch failed
        %             disp(failed.message);
        %             save(sprintf('S04fail-%s.mat',datestr(now,'mmddHHMM')));
        %         end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toBeTried(DD,rossby,JJ,jj)
    [EE,skip]=work_day(DD,JJ(jj),rossby);
    if skip,disp(['skipping ' num2str(jj)]);return;end
    %% save
    save_eddies(EE);
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
    [EE.anticyclones,EE.pass.ac]=walkThroughContsVertically(ee,rossby,cut,DD,-1);
    %% cyclones
    [EE.cyclones,EE.pass.c]=walkThroughContsVertically(ee,rossby,cut,DD,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eddies, pass]=walkThroughContsVertically(ee,rossby,cut,DD,sense)
    pp=0;  pass=initPass(numel(ee))    ;
    %% init
    [eddyType,Zloop]=determineSense(DD.FieldKeys.senses,sense,numel(ee));
    %% loop
    Tv=disp_progress('init','running through conts vertically');
    for kk=Zloop % dir dep. on sense
        Tv= disp_progress('disp',Tv,numel(Zloop),10,1);
        [pass(kk),ee_out]=run_eddy_checks(pass(kk),ee(kk),rossby,cut,DD,sense);
        if all(struct2array(pass(kk))), pp=pp+1;
            %% append healthy found eddy
            eddies(pp)=ee_out;
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
    [ee.circum.si]=EDDyCircumference(zoom);
    %% filter eddies not circle-like enough
    [pass.CR_Shape,ee.isoper, ee.chelt]=CR_Shape(zoom,ee,DD.thresh.shape,DD.switchs);
    if ~pass.CR_Shape, return, end;
    %% get peak position and amplitude w.r.t contour
    [pass.CR_AmpPeak,ee.peak,zoom.ssh_BasePos]=CR_AmpPeak(ee,zoom,DD.thresh.amp);
    if ~pass.CR_AmpPeak, return, end;
    %% get profiles
    [ee.profiles]=EDDyProfiles(ee,zoom.fields);
    %% get radius according to max UV ie min vort
    ee.radius=EDDyRadiusFromUV(ee.peak.z, ee.profiles,zoom.fields);
    %% test
    pass.CR_radius=CR_radius(ee.radius.mean,DD.thresh.radius);
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
% function ee=correctXoverlap(ee,DD)
%     X=DD.map.window.fullsize(2);
%     Y=DD.map.window.size.Y;
%     [ee.coordinates.exact.x]=wrapXidx(ee.coordinates.exact.x,X);
%     [ee.coordinates.int.x]=wrapXidx(ee.coordinates.int.x,X);
%     [ee.centroid.x]=wrapXidx(ee.centroid.x,X);
%     [ee.trackref.x]=wrapXidx(ee.trackref.x,X);
%     [ee.volume.center.x]=wrapXidx(ee.volume.center.x,X);
%     %%
%     ee.centroid.lin=drop_2d_to_1d(ee.centroid.y,ee.centroid.x,Y);
%     ee.trackref.lin=drop_2d_to_1d(ee.trackref.y,ee.trackref.x,Y);
%     ee.volume.center.lin=drop_2d_to_1d(ee.volume.center.y,ee.volume.center.x,Y);
%     %%
%     function [data]=wrapXidx(data,X)
%         data(data<0.5)=X;
%         data(data<1)=1;
%         needcorr=(data>X);
%         data(needcorr)=data(needcorr)-X;
%         data(data<0.5)=X;
%         data(data<1)=1;
%         data(data>X)=X;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks
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
function pass=CR_radius(radius,thresh)
    if radius>=thresh, pass=true; else pass=false; end
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
    [peak.z.y,peak.z.x]=raise_1d_to_2d(diff(z.limits.y)+1, peak.lin);
    peak.amp.to_mean = z.fields.ssh(peak.lin)-peak.mean_ssh;
    %% coordinates in full map
    peak.y=peak.z.y+z.limits.y(1) -1;
    peak.x=peak.z.x+z.limits.x(1) -1;
    if peak.amp.to_contour>=thresh,	pass=true; 	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass,IQ,chelt]=CR_Shape(z,ee,thresh,switches)
    [passes.iq,IQ]=IsopQuo(ee,thresh.iq);
    [passes.chelt,chelt]=chelton_shape(z);
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
function [pass,chelt]=chelton_shape(z)
    % (diameter of circle with equal area)/(maximum distance between nodes)
    %% get max dist in x | y
    x.min=min(z.coor.int.x);
    y.min=min(z.coor.int.y);
    x.max=max(z.coor.int.x);
    y.max=max(z.coor.int.y);
    maxDist=max([sum(z.fields.DX(x.min:x.max)) sum(z.fields.DY(y.min:y.max))]);
    mlat=abs(nanmean(z.fields.lat(:))) ;
    if mlat> 25
        chelt  = 1 - maxDist/4e5;
    else
        chelt  =  1 - maxDist/(8e5*(25 - mlat)/25 + 4e5);
    end
    if chelt >= 0, pass=true; else pass=false; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pass,isoper]=IsopQuo(ee,thresh)
    %% isoperimetric quotient
    % The isoperimetric quotient of a closed curve is defined as the ratio of the curve area to the area of a circle with same perimeter
    % ie isoper=4pi area/circum^2.  isoper(circle)==1;
    isoper=4*pi*ee.area.total/ee.circum.si^2;
    if isoper >= thresh, pass=true; else pass=false; end
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
    area.RadiusOverRossbyL=sqrt(area.total/pi)/rossbyL;
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
function prof=EDDyProfiles(ee,fields)
    %% detect meridional and zonal profiles shifted to baselevel of current level
    offset_term=ee.peak.amp.to_contour*ee.sense.num-ee.level;
    %%	zonal cut
    prof.x.ssh=fields.ssh(ee.peak.z.y,:) + offset_term;
    prof.x.U=fields.U(ee.peak.z.y,:) ;
    prof.x.V=fields.V(ee.peak.z.y,:) ;
    %% meridional cut
    prof.y.ssh=fields.ssh(:,ee.peak.z.x) + offset_term;
    prof.y.U=fields.U(:,ee.peak.z.x) ;
    prof.y.V=fields.V(:,ee.peak.z.x) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radius=EDDyRadiusFromUV(peak,prof,fields)
    %% differentiate velocities to find first local extremata away from peak ie
    %% maximum orbital speed
    %% ie those distances at which the orbital velocity seizes to increase
    
    cb=@(in) [nan;  reshape(in,[],1) ; nan];
    sm=@(in) smooth(in,5,'lowess');
    di=@(in) diff(in,2);
    x.Vdiff=cb(di(sm(prof.x.ssh)));
    y.Udiff=cb(di(sm(prof.y.ssh)));
    %% sign of angular frequency; flowing north, east of peak and south, west of peak
    sense=sign(x.Vdiff(peak.x)) ;
    %% split into both directions
    coor.Xwest=find(x.Vdiff(1:peak.x)*-sense>=0,1,'last');
    if isempty(coor.Xwest), coor.Xwest=1; end
    coor.Xeast=find(x.Vdiff(peak.x:end)*-sense>=0,1,'first') ;
    if isempty(coor.Xeast)
        coor.Xeast=length(x.Vdiff);
    else
        coor.Xeast=coor.Xeast+peak.x-1;
    end
    %%
    coor.Ysouth=find(y.Udiff(1:peak.y)*(-sense)>=0,1,'last');
    
    if isempty(coor.Ysouth), coor.Ysouth=1; end
    coor.Ynorth=find(y.Udiff(peak.y:end)*(-sense) >=0,1,'first');
    if isempty(coor.Ynorth)
        coor.Ynorth=length(y.Udiff);
    else
        coor.Ynorth=coor.Ynorth+peak.y-1;
    end
    %% radius
    radius.zonal=sum(fields.DX(coor.Xwest+1:coor.Xeast))/2;
    radius.meridional=sum(fields.DY(coor.Ysouth+1:coor.Ynorth))/2;
    radius.mean=mean(struct2array(radius));
    %%
    radius.coor=coor;
    
    %     clf
    %nrm=@(in) (in-min(in))/max(in-min(in));
    %     plot(nrm(prof.y.ssh));
    %     hold on
    %     plot(nrm(sm(prof.y.ssh)),'r')
    %     plot(nrm( y.Udiff),'black')
    %     plot(.5*ones(size( y.Udiff)),'y')
    %     grid on
    %     axis tight
    %     legend('raw ssh','lowess filtered','2nd diff')
    %     plot([double(coor.Ysouth)+.5 double(coor.Ynorth)-.5],[.5 .5],'r*')
    %     set(gca,'ytick',.1:.1:.9,'xticklabel',[],'yticklabel',[])
    %      set(gca,'xtick',1:1:49)
    %     savefig('./',100,1200,600,'prof')
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
function [circum]=EDDyCircumference(z)
    %% hypot exact coor diffs times dxdy at those coors
    x=z.coor.exact.x;
    y=z.coor.exact.y;
    xi=z.coor.int.x;
    yi=z.coor.int.y;
    circum=sum(hypot(diff(y).*z.fields.DY(yi(2:end)),diff(x).*z.fields.DX(xi(2:end))));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask=EDDyCut_maskOld(zoom)
    ndgFromLim = @(lim) ndgrid(lim.y(1):lim.y(2),lim.x(1):lim.x(2)) ;
    msk=@(XX,YY,coor) inpolygon(XX,YY,coor.x,coor.y);
    minmax=@(c) [min(c) max(c)];
    %%
    tightLims.x=minmax(zoom.coor.int.x);
    tightLims.y=minmax(zoom.coor.int.y);
    dummymask=false(size(zoom.fields.ssh));
    [YY,XX]=ndgFromLim(tightLims);
    YXlin=drop_2d_to_1d(YY,XX,size(dummymask,1));
    [mfilled, mrim_only] = msk(XX,YY,zoom.coor.int);
    mask.filled  =buildmask(dummymask,YXlin(mfilled));
    mask.rim_only=buildmask(dummymask,YXlin(mrim_only));
    mask.inside= mask.filled & ~mask.rim_only;
    [mask.size.Y, mask.size.X]=size(dummymask);
    function dummy=buildmask(dummy,xlin)
        dummy(xlin)=true;
    end
end
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
function fields_out=EDDyCut_init(fields_in,zoom)
    ya=zoom.limits.y(1);
    yb=zoom.limits.y(2);
    xa=zoom.limits.x(1);
    xb=zoom.limits.x(2);
    %% cut all fields
    for ff=fieldnames(fields_in)'
        field=ff{1};
        fields_out.(field)=fields_in.(field)(ya:yb,xa:xb);
    end
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
    z.limits=enlarge_window(z.limits,enlargeFac,map.sizePlus) ;
    z.coor.int.x=z.coor.int.x-z.limits.x(1) +1;
    z.coor.int.y=z.coor.int.y-z.limits.y(1) +1;
    z.coor.exact.x=z.coor.exact.x -double(z.limits.x(1))  +1;
    z.coor.exact.y=z.coor.exact.y -double(z.limits.y(1))  +1;
    %%
    if strcmp(map.type,'globe')
        %% in global case dismiss eddies touching zonal boundaries (another copy of these eddies exists that is not touching boundaries, due to the zonal appendage in S00b
        pass(1) =  z.limits.x(1)~=1;
        pass(2) =  z.limits.x(2)~=map.sizePlus.X;
        %% also dismiss eddies the zoom window of which is entirely inside the appended stripe ie beyond X_real(2). another legitimate copy of these exists in the western part of the actual map.)
        pass(3) =  ~all(coor.int.x > map.size.X);
    end
    passout=all(pass);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inout=enlarge_window(inout,factor,dim)
    half_width =round((diff(inout.x)+1)*(factor-1)/2);
    half_height=round((diff(inout.y)+1)*(factor-1)/2);
    inout.x(1)=max([1 inout.x(1)-half_width]);
    inout.x(2)=min([dim.X inout.x(2)+half_width]);
    inout.y(1)=max([1 inout.y(1)-half_height]);
    inout.y(2)=min([dim.Y inout.y(2)+half_height]);
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
            EE(cc).coordinates.int.x=int32(round(EE(cc).coordinates.exact.x));
            EE(cc).coordinates.int.y=int32(round(EE(cc).coordinates.exact.y));
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
        %% the following also takes care of the overlap from S00 in the global case
        % x(x>cut.window.size.X)= x(x>cut.window.size.X)-cut.window.size.X ;
        x(x>cut.dim.X)=cut.dim.X;
        y(y>cut.dim.Y)=cut.dim.Y;
        x(x<1)=1;
        y(y<1)=1;
        ee(jj).coordinates.int.x=x;
        ee(jj).coordinates.int.y=y;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plots4debug(zoom,ee) %#ok<*DEFNU>
    %     close all
    %     figure('Position',[1926 250 560 420])
    subplot(1,2,1)
    contourf(zoom.fields.ssh)
    title(['area:',num2str(round(ee.area.total/1000000)),'  -isop:',num2str(round(ee.isoper*1000)/1000)])
    hold on
    plot(zoom.coor.int.x,zoom.coor.int.y)
    subplot(1,2,2)
    pcolor(zoom.fields.ssh)
    hold on
    plot(zoom.coor.exact.x,zoom.coor.exact.y)
    %     print('-dpng','-r100', [num2str(ee.circum.si),'-',num2str(ee.level)])
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
