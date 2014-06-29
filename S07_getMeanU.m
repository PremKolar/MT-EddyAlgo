%%%%%%%%%
% Created: 08-Apr-2014 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S07_getMeanU
	%% init
	DD=initialise;
	%% find files
	[file]=findVelFiles(DD);
	%% get dims
	[d,pos,dim]=getDims(file,DD);
	%% means
	means=getMeans(d,pos,dim,file,DD); %#ok<NASGU>
	%% save
	save([DD.path.meanU.file], 'means')
	disp(['done!'])
	conclude(DD);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function def=deformation(g)
	%% calc U gradients
	[Z,Y,X]=size(g.u);
	DY=shiftdim(repmat(g.dy,[1,1,Z]),2);
	DX=shiftdim(repmat(g.dx,[1,1,Z]),2);
	def.dudy=cat(2, nan(Z,1,X), diff(g.u,1,2)) ./ DY;
	def.dvdx=cat(3, nan(Z,Y,1), diff(g.v,1,3)) ./ DX;
	def.dvdy=cat(2, nan(Z,1,X), diff(g.v,1,2)) ./ DY;
	def.dudx=cat(3, nan(Z,Y,1), diff(g.u,1,3)) ./ DX;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ow]=getOW(g)
	ow.vorticity = g.def.dvdx - g.def.dudy;
	ow.divergence= g.def.dudx + g.def.dvdy;
	ow.stretch   = g.def.dudx - g.def.dvdy;
	ow.shear     = g.def.dvdx + g.def.dudy;
	%% okubo weiss
	ow.ow=.5*(-ow.vorticity.*2+ow.divergence.*2+ow.stretch.*2+ow.shear.*2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grids=readGrids(file,DD,dim)
	disp(['found ' file.U ' and ' file.V])
	u=squeeze(nc_varget(file.U,DD.map.in.keys.U,dim.start,dim.length))/DD.parameters.meanUunit;
	v=squeeze(nc_varget(file.V,DD.map.in.keys.V,dim.start,dim.length))/DD.parameters.meanUunit;
	[Z,~,~]=size(u);
	for z=1:Z
		grids.u(z,:,:)=squeeze(u(z,:,:))-smooth2a(squeeze(u(z,:,:)),20);
		grids.v(z,:,:)=squeeze(v(z,:,:))-smooth2a(squeeze(v(z,:,:)),20);
	end	
	grids.lat=nc_varget(file.V,DD.map.in.keys.lat,dim.start(3:4),dim.length(3:4));
	grids.lon=nc_varget(file.V,DD.map.in.keys.lon,dim.start(3:4),dim.length(3:4));
% 	lat=shiftdim(repmat(lat,[1,1,size(grids.u,1)]),2);
% 	lon=shiftdim(repmat(lon,[1,1,size(grids.u,1)]),2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dy,dx]=dydx(g)
	%% grid increment sizes
	dy=deg2rad(abs(diff(double(g.lat),1,1)))*earthRadius;
	dx=deg2rad(abs(diff(double(g.lon),1,2)))*earthRadius.*cosd(g.lat(:,1:end-1));
	%% append one line/row to have identical size as other fields
	dy=dy([1:end end],:);
	dx=dx(:,[1:end end]);
	%% correct 360° crossings
	seamcrossflag=dx>100*median(dx(:));
	dx(seamcrossflag)=abs(dx(seamcrossflag) - 2*pi*earthRadius.*cosd(g.lat(seamcrossflag)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [means]=GMzFromOWcase(file,DD,dim)
	for kk=1:numel(file)
		%%
		grids=readGrids(file(kk),DD,dim)
% 		%%
% 		fltr=20;
% 		smoothGrids(grids,fltr)
		%%
		[grids.dy,grids.dx]=dydx(grids)
		%% deformation
		grids.def=deformation(grids);
		%%
		ow=getOW(grids);
		%%
		dsfgb
	end
	means=0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function means=getMeans(d,pos,dim,file,DD)
	if DD.switchs.meanUviaOW
		[means]=GMzFromOWcase(file,DD,dim);
	else
		[means]=GMzConstCase(file,DD,dim);
	end
	
	
	
	
	
	
	
	
	
	
	
	%-----------------------------------------------------------------------
	function [means]=GMzConstCase(file,DD,dim)
		for kk=1:numel(file)
			disp(['found ' file(kk).U ' and ' file(kk).V])
			U(:,:,kk)=squeeze(nc_varget(file(kk).U,DD.map.in.keys.U,dim.start,dim.length))/DD.parameters.meanUunit; %#ok<*AGROW>
			V(:,:,kk)=squeeze(nc_varget(file(kk).V,DD.map.in.keys.V,dim.start,dim.length))/DD.parameters.meanUunit;
			%%
			x=DD.map.window.size.X;
			y=DD.map.window.size.Y;
			if x~=size(U,2) || y~=size(U,1)
				warning('trivially resizing U/V data!!! ') %#ok<WNTAG>
				sleep(5)
				U=downsize(U,x,y);
				V=downsize(V,x,y);
			end
		end
		disp(['creating means'])
		U(U<-1e33)=nan; % missing values
		V(V<-1e33)=nan; % missing values
		means.zonal=nanmean(U,3);
		means.merid=nanmean(V,3);
		means.total=hypot(means.zonal,means.merid);
		means.direc=azimuth(zeros(size(means.zonal)),zeros(size(means.zonal)),means.merid,means.zonal);
		means.depth=d(pos.z.start);
		%%
		disp(['resizing to output size'])
		lin=extractfield(load(DD.path.protoMaps.file),'idx');
		means.small.zonal=nan(DD.map.out.Y,DD.map.out.X);
		for li=unique(lin(lin~=0 & ~isnan(lin)))
			means.small.zonal(li)=nanmean(means.zonal(lin==li));
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,pos,dim]=getDims(file,DD)
	if DD.switchs.meanUviaOW
		d=nan;
		[pos,dim]=DIMzFromOWcase(DD);
	else
		[d,pos,dim]=DIMzConstCase(file,DD);
	end
	%-----------------------------------------------------------------------
	function [d,pos,dim]=DIMzConstCase(file,DD)
		dWanted=DD.parameters.meanU;
		d=nc_varget(file(1).U,DD.map.in.keys.z);
		[~,pos.z.start]=min(abs(d-dWanted));
		pos.z.start=pos.z.start - 1; % starts at 0
		pos.z.length=1;
		pos.x.start=DD.map.window.limits.west - 1;
		pos.x.length=DD.map.window.size.X;
		pos.y.start=DD.map.window.limits.south-1;
		pos.y.length=DD.map.window.size.Y;
		dim.start = [0 pos.z.start pos.y.start pos.x.start];
		dim.length = 	[inf pos.z.length pos.y.length pos.x.length ];
	end
	%-----------------------------------------------------------------------
	function [pos,dim]=DIMzFromOWcase(DD)
		pos.z.start=0;
		pos.z.length=inf;
		pos.x.start=DD.map.window.limits.west - 1;
		pos.x.length=DD.map.window.size.X;
		pos.y.start=DD.map.window.limits.south-1;
		pos.y.length=DD.map.window.size.Y;
		dim.start = [0 pos.z.start pos.y.start pos.x.start];
		dim.length = 	[inf pos.z.length pos.y.length pos.x.length ];
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [file]=findVelFiles(DD)
	%% find the U and V files
	ucc=0; vcc=0;
	file=struct;
	for kk=1:numel(DD.path.raw.files)
		if ~isempty(strfind(DD.path.raw.files(kk).name,'UVEL'))
			ucc=ucc+1;
			file(ucc).U=[DD.path.raw.name DD.path.raw.files(kk).name]; %#ok<AGROW>
		end
		if ~isempty(strfind(DD.path.raw.files(kk).name,'VVEL'))
			vcc=vcc+1;
			file(vcc).V=[DD.path.raw.name DD.path.raw.files(kk).name]; %#ok<AGROW>
		end
	end
	if isempty(file)
		disp(['put U/V files into ' DD.path.raw])
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
