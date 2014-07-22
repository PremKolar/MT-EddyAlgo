%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 16-Jul-2014 13:52:44
% Computer:GLNX86
% Matlab:7.9
% Author:NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOWmain(DD)
	[dim,raw] = setup(DD);
	if DD.debugmode
		spmd_body(DD,dim,raw);
	else
		spmd(DD.threads.num)
			spmd_body(DD,dim,raw);
			% 			disp_progress('conclude');
		end
	end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD,dim,raw)
	id = labindex;
	lims = DD.TSow.lims-1;
	timesteps = lims(id,1):lims(id,2);
	keys = DD.TS.keys;
	fn=0;
	Ysteps=round(linspace(0,dim.ws.Y-1,10))
	for cc =1:numel(Ysteps)-1
		Ychunks{cc}=Ysteps(cc):Ysteps(cc+1)-1;
	end
	
	
	
	for tt=timesteps
		fn=fn+1;
		fname(fn)={initNcFile('density',dim.ws,tt,DD.path.TSow)};
	end
	fnameRhoMeanPart=initNcRhoSumPart('densitySum',dim.ws,id)
	
	
	
	
	%%
	
	for tt = timesteps
		Fin = DD.path.TSow.files(tt+1);
		for cc = numel(Ysteps)-1
			yy=Ysteps(cc);
			Ychunk=Ychunks{cc};
			dim.len(2) = length(Ychunk);
			dim.strt(2) = yy;
			[TS] = getNowAtY(Fin,keys,dim,DD.path.full3d.name,raw,Ychunks{cc}+1);
			rho = calcDens(TS,dim.len);
			nc_varput(fname{tt+1},'density',rho,dim.out.strt,dim.len);
		end
	end
	
	
	%%
	for tt = timesteps
		Fin = DD.path.TSow.files(tt+1);
		for cc = 1:numel(Ysteps)-1
			yy=Ysteps(cc);
			Ychunk=Ychunks{cc};
			dim.len(2) = length(Ychunk);
			dim.strt(2) = yy;
rho=nc_varget(fname{tt+1},'density',rho,dim.out.strt,dim.len)		




	[TS] = getNowAtY(Fin,keys,dim,DD.path.full3d.name,raw,Ychunks{cc}+1);
			rho = calcDens(TS,dim.len);
			nc_varput(fname{tt+1},'density',rho,dim.out.strt,dim.len);
		end
	end
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = getNowAtY(Fin,keys,dim,dirOut,raw,YY)
	make1d = @(Axb,C) reshape(repmat(double(Axb),C,1),[],1)   ;
	Z=dim.ws.Z;
	out.lat = make1d(raw.lat(YY,:),Z);
	out.lon = make1d(raw.lon(YY,:),Z);
	out.dx = make1d(raw.dx(YY,:),Z);
	out.dy = make1d(raw.dy(YY,:),Z);
	out.depth = repmat(double(raw.depth),numel(YY),1);
	out.temp = double(sparse(reshape(nc_varget(Fin.temp,keys.temp,dim.strt, dim.len),[],1)));
	out.salt = double(sparse(reshape(nc_varget(Fin.salt,keys.salt,dim.strt, dim.len),[],1)))*1000;
	out.file = fileStuff(Fin.salt, dirOut);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RHO] = calcDens(TS,dim)
	RHO=reshape(dens(TS),dim)   ;
	function rho=dens(TS)
		rho=sw_dens(TS.salt,TS.temp,sw_pres(TS.depth,TS.lat));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [dim,raw]=setup(DD)
	ws = DD.TSow.window.size;
	wl = DD.TSow.window.limits;
	Fin = DD.path.TSow.files(1);
	keys = DD.TS.keys;
	dimConst.strt = [ wl.south-1 wl.west-1 ];
	dimConst.len = [ ws.Y ws.X ];
	raw.depth = nc_varget(Fin.salt,keys.depth);
	raw.lat = nc_varget(Fin.temp, keys.lat, dimConst.strt, dimConst.len);
	raw.lon = nc_varget(Fin.temp, keys.lon, dimConst.strt, dimConst.len);
	[raw.dy,raw.dx] = getdydx( raw.lat, raw.lon);
	raw.corio = coriolisStuff(raw.lat);
	dim.strt(1) = 0;
	dim.strt(2:3) = dimConst.strt;
	dim.len(1) = ws.Z;
	dim.len(3) = dimConst.len(2);
	%     dim.len(2) = 1;
	dim.ws=ws;
	dim.out.strt=[0 0 0];
	dim.out.len=[1 dimConst.len];
end

function meanDensInTime(raw)
	DD.path.full3d.name
	
	
	if exist(file_out,'file');
		dispM('exists');
		%return;
	end
	
	
	%% OW
	[raw.OW] = calcOW(raw,cc);
	%% save
	dispM('saving..')
	saveChunk(raw,file_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Calculations(DD,ff,raw)
	
	%% merge
	file_out = [DD.path.Rossby.name,'OW_',sprintf('%03d',ff),'.mat'];
	if exist(file_out,'file');
		dispM('exists');
		%return;
	end
	raw = getNowAtX(raw,DD,ff);
	
	%% calculate Densisty
	[raw.rho] = calcDens(raw,cc);
	%% OW
	[raw.OW] = calcOW(raw,cc);
	%% save
	dispM('saving..')
	saveChunk(raw,file_out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OW = 	calcOW(raw,cc)
	dispM(['getting okubo weiss ',cc],1)
	%% rho gradient
	[gr.drdx,gr.drdy] = getDrhodx(raw.rho,raw.dx,raw.dy);
	%% velocities
	vels = getVels(raw,gr);
	clear gr;
	%% uvgrads
	uvg = UVgrads(vels,raw.dx,raw.dy);
	clear vels;
	%% deformation
	deform = getDefo(uvg);
	clear uvg;
	%% okubo weiss
	OW = okuweiss(deform);
end
%% ---------------------------------------------------------------------
function OW = okuweiss(d)
	OW = (-d.vorticity.*2+d.divergence.*2+d.stretch.*2+d.shear.*2)/2;
	OW(abs(OW)>1) = nan;
end
%-----------------------------------------------------------------------
function defo = getDefo(uvg)
	defo.vorticity = uvg.dVdx - uvg.dUdy;
	defo.divergence = uvg.dUdx + uvg.dVdy;
	defo.stretch = uvg.dUdx - uvg.dVdy;
	defo.shear = uvg.dVdx + uvg.dUdy;
end
%-----------------------------------------------------------------------
function uvg = UVgrads(vels,dx,dy)
	%% calc U gradients
	dUdy = diff(vels.U,1,2);
	dUdx = diff(vels.U,1,3);
	dVdy = diff(vels.V,1,2);
	dVdx = diff(vels.V,1,3);
	[Z,~,~] = size(vels.U);
	uvg.dUdy = dUdy(:, [1:end, end], :)./ vertstack(dy,Z);
	uvg.dUdx = dUdx(:, :,[1:end, end] )./ vertstack(dx,Z);
	uvg.dVdy = dVdy(:, [1:end, end], :)./ vertstack(dy,Z);
	uvg.dVdx = dVdx(:, :,[1:end, end] )./ vertstack(dx,Z);
end
%---------
function vels = getVels(raw,gr)
	cor = raw.corio;
	depth = raw.depth;
	rhoRef = 1000;
	[Z,Y,X] = size(gr.drdy);
	gzOverRhoF = vertstack(cor.GOverF,Z) .* repmat(depth,[1,Y,X]) / rhoRef;
	vels.U = -gr.drdy .* gzOverRhoF;
	vels.V = gr.drdx .* gzOverRhoF;
	% semi.x = raw.filterRadius;
	% semi.y = raw.filterRadius;
	% T = disp_progress('init','high pass filtering geostrophic velocities');
	% for z = 1:Z
	% T = disp_progress('init',T,Z,10);
	% [vels.U(z,:,:)] = ellipseFltr(semi,squeeze(U(z,:,:)));
	% [vels.V(z,:,:)] = ellipseFltr(semi,squeeze(V(z,:,:)));
	% end
	
	% for z = 1:Z
	% vels.U(z,:,:) = smooth2a(squeeze(U(z,:,:)),3);
	% vels.V(z,:,:) = smooth2a(squeeze(V(z,:,:)),3);
	% end
	% for y = 1:Y
	% Uy(:,y,:) = smooth2a(squeeze(U(:,y,:)),3);
	% Vy(:,y,:) = smooth2a(squeeze(V(:,y,:)),3);
	% end
	% for x = 1:X
	% Ux(:,:,x) = smooth2a(squeeze(U(:,:,x)),3);
	% Vx(:,:,x) = smooth2a(squeeze(V(:,:,x)),3);
	% end
	% vels.U = reshape(mean([Uz(:),Uy(:),Ux(:)],2),[Z,Y,X]);
	% vels.V = reshape(mean([Uz(:),Uy(:),Ux(:)],2),[Z,Y,X]);
	
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [drdx,drdy] = getDrhodx(rho,dx,dy)
	%% calc density gradients
	[Z,~,~] = size(rho);
	drdx = diff(rho,1,3);
	drdy = diff(rho,1,2);
	drdx = drdx(:,:,[1:end, end]) ./ vertstack(dx,Z);
	drdy = drdy(:,[1:end, end],:) ./ vertstack(dy,Z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = fileStuff(fin,dirout)
	dateIdx = regexp(fin,'[0-9]{8}');
	dateN = fin(dateIdx:dateIdx+7);
	out.date = datenum(dateN,'yyyymmdd');
	out.fout = [ dirout, 'meanDens',dateN,'.mat'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = coriolisStuff(lat)
	OmegaTw = 2*angularFreqEarth;
	%% f
	out.f = OmegaTw*sind(lat);
	%% beta
	out.beta = OmegaTw/earthRadius*cosd(lat);
	%% gravity
	g = sw_g(lat,zeros(size(lat)));
	%% g/f
	out.GOverF = g./out.f;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dy,dx] = getdydx(lat,lon)
	%% grid increment sizes
	dy = deg2rad(abs(diff(double(lat),1,1)))*earthRadius;
	dx = deg2rad(abs(diff(double(lon),1,2)))*earthRadius.*cosd(lat(:,1:end-1));
	%% append one line/row to have identical size as other fields
	dy = dy([1:end end],:);
	dx = dx(:,[1:end end]);
	%% correct 360° crossings
	seamcrossflag = dx>100*median(dx(:));
	dx(seamcrossflag) = abs(dx(seamcrossflag) - 2*pi*earthRadius.*cosd(lat(seamcrossflag)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveChunk(raw,file_out) %#ok<INUSL>
	save(file_out,'-struct','raw');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function name=initNcRhoSumPart(toAdd,WinSize,id)
	name.sum=['./threadChunkRhoSum' sprintf('%02d.nc',id) ];
	name.count=['./threadChunkRhoCount' sprintf('%02d.nc',id) ];
	nc_create_empty(name.sum,'clobber');
	nc_create_empty(name.count,'clobber');
	nc_adddim(name.sum,'k_index',WinSize.Z);
	nc_adddim(name.sum,'i_index',WinSize.X);
	nc_adddim(name.sum,'j_index',WinSize.Y);
	nc_adddim(name.count,'k_index',WinSize.Z);
	nc_adddim(name.count,'i_index',WinSize.X);
	nc_adddim(name.count,'j_index',WinSize.Y);
	%% rho summed
	varstruct.Name = toAdd;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'k_index','j_index','i_index' };
	nc_addvar(name.sum,varstruct)
	%% count
	varstruct.Name = [toAdd '_count'];
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'k_index','j_index','i_index' };
	nc_addvar(name.count,varstruct)	
end

function fname=initNcFile(toAdd,WinSize,ff,ds)
	fname=[ds.dailyBaseName sprintf('%04d.nc',ff) ];
	nc_create_empty(fname,'clobber');
	nc_adddim(fname,'k_index',WinSize.Z);
	nc_adddim(fname,'i_index',WinSize.X);
	nc_adddim(fname,'j_index',WinSize.Y);
	%% rho
	varstruct.Name = toAdd;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'k_index','j_index','i_index' };
	nc_addvar(fname,varstruct)
	%% depth
	varstruct.Name = 'depth';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'k_index'};
	nc_addvar(fname,varstruct)
	%% lat
	varstruct.Name = 'lat';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'j_index','i_index' };
	nc_addvar(fname,varstruct)
	%% lon
	varstruct.Name = 'lon';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'j_index','i_index' };
	nc_addvar(fname,varstruct)
end