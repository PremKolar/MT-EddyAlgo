%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% walks through all the contours and decides whether they qualify 
function S03_filter_eddies
	%% init
	DD=initialise('conts');
	init_threads(DD.threads.num);
	%% spmd
	%spmd
		id=labindex;
		spmd_body(DD,id)
	%end
	%% update infofile
	save_info(DD)
end
%% main functions
function save_eddies(EE)
	save(EE.filename,'-struct','EE')
end
function spmd_body(DD,id)
	Td=disp_progress('init','filtering contours');
	for jj=DD.threads.lims(id,1):DD.threads.lims(id,2)
		Td=disp_progress('disp',Td,diff(DD.threads.lims(id,:))+1,4242);
		%%
		EE=work_day(DD,jj);
		%% save
		save_eddies(EE);
	end
end
function EE=work_day(DD,jj)
	%% get ssh data
	cut=read_fields(DD,jj,'cuts');
	%% get contours
	cont=read_fields(DD,jj,'conts');
	%% put all eddies into a struct: ee(number of eddies).characteristica  
	ee=eddies2struct(cont.all,DD.thresh.corners);
	%% avoid out of bounds integer coordinates close to boundaries
	[ee_clean,cut]=CleanEddies(ee,cut);
	%% find them
	EE=find_eddies(ee_clean,cut,DD,jj);
end
function EE=find_eddies(ee_clean,cut,DD,jj)
	%% init
	EE=struct;
	EE.filename=[DD.path.eddies.name, regexprep(DD.path.conts.files(jj).name, 'CONT', 'EDDIE')];
	%% anti cyclones
	EE.anticyclones=anti_cyclones(ee_clean,cut,DD);
	%% cyclones
	EE.cyclones=cyclones(ee_clean,cut,DD);
end
function ACyc=anti_cyclones(ee,cut,DD)
	PASS=false(numel(ee),1);	pp=0;
	%% loop over eddies, starting at deepest eddies, upwards
	for kk=1:numel(ee)
		[PASS(kk),ee_out]=run_eddy_checks(ee(kk),cut,DD,-1);
		if PASS(kk), pp=pp+1;
			%% append healthy found eddy
			ACyc(pp)=ee_out;  %#ok<AGROW>
			%% nan out ssh where eddy was found
			cut.grids.SSH(ee_out.mask)=nan;						
		end
	end
	
	if ~any(PASS)
		error('no anticyclones made it through the filter...')
	end
end
function Cyc=cyclones(ee,cut,DD)
	PASS=false(numel(ee),1);	pp=0;
	%% loop over eddies, starting at highest eddies, downwards
	for kk=numel(ee):-1:1
		[PASS(kk),ee_out]=run_eddy_checks(ee(kk),cut,DD,1);
		if PASS(kk),	pp=pp+1;
			%% append healthy found eddy
			Cyc(pp)=ee_out;
			%% nan out ssh where eddy was found
			cut.grids.SSH(ee_out.mask)=nan;
		end
	end
	if ~any(PASS)
		error('no cyclones made it through the filter...')
	end
end
function [pass,ee]=run_eddy_checks(ee,cut,dd,direction)
	%% pre-nan-check
	pass=CR_RimNan(ee.coordinates.int, cut.dim.Y	, cut.grids.SSH);
	if ~pass, return, end;
	%% corners-check
	pass=CR_corners(ee.circum.length	,dd.thresh.corners);
	if ~pass, return, end;
	%% closed ring check
	[pass]=CR_ClosedRing(ee);
	if ~pass, return, end;
	%% pre filter 'thin 1dimensional' eddies
	pass=CR_2dEddy(ee.coordinates.int);
	if ~pass, return, end;
	%% get coordinates for zoom cut
	[zoom]=get_window_limits(ee.coordinates,cut.dim,4);
	%% cut out rectangle encompassing eddy range only for further calcs
	zoom.fields=EddyCut_init(cut.grids,zoom);  
	%% generate logical masks defining eddy interiour and outline
	zoom.mask=EddyCut_mask(zoom);
	%% check for nans within eddy
	[pass]=CR_Nan(zoom);
	if ~pass, return, end;
	%% check for correct sense
	[pass,ee.sense]=CR_sense(zoom,direction,ee.level);
	if ~pass, return, end;
	%% calculate area with respect to contour
	[ee.area]=Area(zoom);
	%% calc contour circumference in [SI]
	[ee.circum.si]=EddyCircumference(zoom);
	%% filter eddies not circle-like enough
	[pass,ee.isoper]=CR_Shape(ee,dd.thresh.shape,dd.switchs.IQ);
	if ~pass, return, end;
	%% get peak position and amplitude w.r.t contour
	[pass,ee.peak,zoom.SSH_BasePos]=CR_AmpPeak(ee,zoom,dd.thresh.amp);
	if ~pass, return, end;
	%% get profiles
	[ee.profiles]=EddyProfiles(ee,zoom.fields);
	%% get radius according to max UV ie min vort
	ee.radius=EddyRadiusFromUV(ee.peak.z, ee.profiles,zoom.fields);
	%% test
	pass=CR_radius(ee.radius.mean,dd.thresh.radius);
	if ~pass, return, end;
	%% get ideal ellipse contour
	zoom.mask.ellipse=EddyEllipse(ee,zoom.mask);
	%% get effective amplitude relative to ellipse;
	[pass,ee.peak.amp.to_ellipse]=EddyAmp2Ellipse(ee.peak.lin,zoom,dd.thresh.amp);
	if ~pass, return, end;
	%% get area regarding ellipse;
	%redundant!
	%[ee.area.ellipse]=EddyArea2Ellipse(ee.radius.mean);
	%% append mask to ee in cut coordinates
	[ee.mask]=sparse(EddyPackMask(zoom.mask.filled,zoom.limits,size(cut.grids.SSH)));
	%plots4debug(zoom,ee)
	%% get center of 'volume'
	[ee.volume]=CenterOfVolume(zoom,ee.area.total,cut.dim.Y);
	%% get coordinates
	[ee.geo]=geocoor(zoom,ee.volume);
	%% append 'age'
	ee.age=0;	
end
%% checks
function [pass,sense]=CR_sense(zoom,direc,level)
	pass=false;
	sense=struct;
	%% water column up: seeking anti cyclones; down: cyclones
	if direc==-1
		if all(zoom.fields.SSH(zoom.mask.inside) >= level )
			pass=true;
			sense.str='AntiCyclonic';
			sense.num=-1;
		end
	elseif direc==1
		if all(zoom.fields.SSH(zoom.mask.inside) <= level )
			pass=true;
			sense.str='Cyclonic';
			sense.num=1;
		end
	end
end
function pass=CR_radius(radius,thresh)
	if radius>=thresh, pass=true; else pass=false; end
end
function pass=CR_RimNan(coor, Y, SSH)
	pass=true;
	if any(isnan(SSH(drop_2d_to_1d(coor.y, coor.x, Y)))), pass=false; end
end
function pass=CR_corners(corners,thresh)
	pass=true;
	if corners < thresh, pass=false; end
end

function [pass,peak,base]=CR_AmpPeak(ee,z,thresh)
	pass=false;
	%%
	peak.mean_SSH=mean(z.fields.SSH(z.mask.filled));
	%% make current level zero level and zero out everything else
	base=poslin(-ee.sense.num*(z.fields.SSH-ee.level));
	base(~z.mask.filled)=0;
	%% amplitude
	[peak.amp.to_contour,peak.lin]=max(base(:));
	[peak.z.y,peak.z.x]=raise_1d_to_2d(diff(z.limits.y)+1, peak.lin);
	peak.amp.to_mean.of_contour=peak.amp.to_contour - mean(base(:));
	%% coordinates in full map
	peak.y=peak.z.y+z.limits.y -1;
	peak.x=peak.z.x+z.limits.x -1;
	if peak.amp.to_contour>=thresh,	pass=true; 	end
end
function [pass,ShapeRatio]=CR_Shape(ee,thresh,iq)
if iq 	
	[pass,ShapeRatio]=IsopQuo(ee,thresh.iq);
else
	[pass,ShapeRatio]=chelton_shape(ee,thresh.chelt);
end
end

function [pass,chelt]=chelton_shape(ee,thresh)
% (diameter of circle with equal area)/(maximum distance between nodes)
%% get max dist in x | y
ihgkjh
if chelt >= thresh, pass=true; else pass=false; end
end


function [pass,isoper]=IsopQuo(ee,thresh)
	%% isoperimetric quotient
% The isoperimetric quotient of a closed curve is defined as the ratio of the curve area to the area of a circle with same perimeter 
% ie isoper=4pi area/circum^2.  isoper(circle)==1; 
	isoper=12.5664*ee.area.total/ee.circum.si^2;
	if isoper >= thresh, pass=true; else pass=false; end
end


function [pass]=CR_2dEddy(coor)
	if (max(coor.x)-min(coor.x)<2) || (max(coor.y)-min(coor.y)<2)
		pass=false;
	else
		pass=true;
	end
end
function [pass]=CR_Nan(z)
	ssh=z.fields.SSH(z.mask.filled);
	if ~any(isnan(ssh(:))), pass=true; else pass=false; end
end
function [pass]=CR_ClosedRing(ee)
	x=ee.coordinates.int.x;
	y=ee.coordinates.int.y;
	if abs(x(1)-x(end))>1 || abs(y(1)-y(end))>1;
		pass=false;
	else
		pass=true;
	end
end





%% others
function [area]=Area(z)
	area=struct;
	area.pixels=(z.fields.DX.*z.fields.DY).*(z.mask.inside + z.mask.rim_only/2);  % include 'half of rim'
	area.total=sum(area.pixels(:));
end
function	[mask_out]=EddyPackMask(mask_in,limits,dims)
	mask_out=false(dims);
	mask_out(limits.y(1):limits.y(2),limits.x(1):limits.x(2))=mask_in;
end


function [pass,amp] = EddyAmp2Ellipse(peak,zoom,thresh)
	%% mean amplitude with respect to ellipse contour 
	amp=abs(zoom.fields.SSH(peak)-nanmean(zoom.fields.SSH(zoom.mask.ellipse)));
	if amp>=thresh, pass=true; else pass=false; end
end
function [ellipse]=EddyEllipse(ee,mask)
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
function prof=EddyProfiles(ee,fields)
	%% detect meridional and zonal profiles shifted to baselevel of current level
	offset_term=ee.peak.amp.to_contour*ee.sense.num-ee.level;
	%%	zonal cut
	prof.x.ssh=fields.SSH(ee.peak.z.y,:) + offset_term;
	prof.x.U=fields.U(ee.peak.z.y,:) ;
	prof.x.V=fields.V(ee.peak.z.y,:) ;
	%% meridional cut
	prof.y.ssh=fields.SSH(:,ee.peak.z.x) + offset_term;
	prof.y.U=fields.U(:,ee.peak.z.x) ;
	prof.y.V=fields.V(:,ee.peak.z.x) ;
end
function radius=EddyRadiusFromUV(peak,prof,fields)
	%% differentiate velocities to find first local extremata away from peak ie
	%% maximum orbital speed
	%% ie those distances at which the orbital velocity seizes to increase 
	x.Vdiff=(diff(smooth(prof.x.V)));
	y.Udiff=(diff(smooth(prof.y.U)));
	%% append one value in case eddy close to wall
	x.Vdiff=[x.Vdiff; x.Vdiff(end)];
	y.Udiff=[y.Udiff; y.Udiff(end)];
	%% sign of angular frequency; flowing north, east of peak and south, west of peak
	sense=sign(x.Vdiff(peak.x)) ;
	%% split into both directions
	coor.Xwest=find(x.Vdiff(1:peak.x)*sense<=0,1,'last');
	if isempty(coor.Xwest), coor.Xwest=1; end
	coor.Xeast=find(x.Vdiff(peak.x:end)*sense<=0,1,'first') ;
	if isempty(coor.Xeast)
		coor.Xeast=length(x.Vdiff);
	else
		coor.Xeast=coor.Xeast+peak.x;
	end
	%%
	coor.Ysouth=find(y.Udiff(1:peak.y)*(-sense)<=0,1,'last');
	if isempty(coor.Ysouth), coor.Ysouth=1; end
	coor.Ynorth=find(y.Udiff(peak.y:end)*(-sense) <=0,1,'first');
	if isempty(coor.Ynorth)
		coor.Ynorth=length(y.Udiff);
	else
		coor.Ynorth=coor.Ynorth+peak.y;
	end
	%% radius
	radius.zonal=sum(fields.DX(coor.Xwest:coor.Xeast))/2;
	radius.meridional=sum(fields.DY(coor.Ysouth:coor.Ynorth))/2;
	radius.mean=mean(struct2array(radius));
	%%
	radius.coor=coor;
end
function [geo]=geocoor(zoom,volume)
	xz=volume.center.xz;
	yz=volume.center.yz;
	geo.lat=interp2(zoom.fields.LAT,xz,yz);
	geo.lon=interp2(zoom.fields.LON,xz,yz);

end
function [volume]=CenterOfVolume(zoom,area,Y)
	%% get "volume" of eddy
	ssh=zoom.SSH_BasePos;
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
end
function [circum]=EddyCircumference(z)
	%% get perimeter coor and check where x or y change
	x=z.coor.int.x	;
	y=z.coor.int.y;
	dx=[0 ;abs(double(logical(diff(x))))];
	dy=[0 ;abs(double(logical(diff(y))))];
	circum=sum(hypot(z.fields.DX(x).*dx,z.fields.DY(y).*dy)); %positive bias at bad resolution!
end
function mask=EddyCut_mask(zoom)
	[Y,X]=size(zoom.fields.SSH);
	mask.rim_only=false(Y,X);
	mask.rim_only(sub2ind([Y,X], zoom.coor.int.y, zoom.coor.int.x))=true;
	xma=max(zoom.coor.int.x);
	xmi=min(zoom.coor.int.x);
	yma=max(zoom.coor.int.y);
	ymi=min(zoom.coor.int.y);
	mask.filled=mask.rim_only;
	mask.filled(ymi:yma,xmi:xma)=logical(imfill(mask.rim_only(ymi:yma,xmi:xma),'holes'));
	mask.inside= mask.filled & ~mask.rim_only;
	mask.size.Y=Y;
	mask.size.X=X;
end
function fields_out=EddyCut_init(fields_in,zoom)
	ya=zoom.limits.y(1);
	yb=zoom.limits.y(2);
	xa=zoom.limits.x(1);
	xb=zoom.limits.x(2);
	for ff=fieldnames(fields_in)'
		field=ff{1};
		fields_out.(field)=fields_in.(field)(ya:yb,xa:xb);
	end
end
function [z]=get_window_limits(coor,dim,enlargeFac)
	z.coor=coor;
	%% output
	z.limits.x(1)=min(coor.int.x);
	z.limits.y(1)=min(coor.int.y);
	z.limits.x(2)=max(coor.int.x);
	z.limits.y(2)=max(coor.int.y);
	z.limits=enlarge_window(z.limits,enlargeFac,dim) ;
	z.coor.int.x=z.coor.int.x-z.limits.x(1) +1;
	z.coor.int.y=z.coor.int.y-z.limits.y(1) +1;
	z.coor.exact.x=z.coor.exact.x -double(z.limits.x(1))  +1;
	z.coor.exact.y=z.coor.exact.y -double(z.limits.y(1))  +1;
end
function inout=enlarge_window(inout,factor,dim)
	half_width =round((diff(inout.x)+1)*(factor-1)/2);
	half_height=round((diff(inout.y)+1)*(factor-1)/2);
	inout.x(1)=max([1 inout.x(1)-half_width]);
	inout.x(2)=min([dim.X inout.x(2)+half_width]);
	inout.y(1)=max([1 inout.y(1)-half_height]);
	inout.y(2)=min([dim.Y inout.y(2)+half_height]);
end
function [EE]=eddies2struct(CC,thresh)
	EE=struct;
	ii=1;cc=0;
	while ii<size(CC,1);
		len=  CC(ii,2);% contourc saves the length of each contour before appending the next
		if len>thresh
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
function [ee,cut]=CleanEddies(ee,cut)
	[cut.dim.Y,cut.dim.X]=size(cut.grids.SSH);
	for jj=1:numel(ee)
		x=ee(jj).coordinates.int.x;
		y=ee(jj).coordinates.int.y;
		x(x<1)=1;
		x(x>cut.dim.X)=cut.dim.X;
		y(y<1)=1;
		y(y>cut.dim.Y)=cut.dim.Y;
		ee(jj).coordinates.int.x=x;
		ee(jj).coordinates.int.y=y;
	end
end
function plots4debug(zoom,ee) %#ok<DEFNU>
	close all
	figure('Position',[1926 250 560 420])
	subplot(1,2,1)
	contourf(zoom.fields.SSH)
	title(['area:',num2str(round(ee.area.total/1000000)),'  -isop:',num2str(round(ee.isoper*1000)/1000)])
	hold on
	plot(zoom.coor.int.x,zoom.coor.int.y)
	subplot(1,2,2)
	pcolor(zoom.fields.SSH)
	hold on
	plot(zoom.coor.exact.x,zoom.coor.exact.y)
	print('-dpng','-r100', [num2str(ee.circum.si),'-',num2str(ee.level)])
end
