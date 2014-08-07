%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 16-Jul-2014 13:52:44
% Computer:GLNX86
% Matlab:7.9
% Author:NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metaD=maxOWmain(DD)
	[Dim,raw] = setup(DD);
	%%
	spmd(DD.threads.num)
		daily=initbuildRho(Dim,DD);
	end
	buildRho(daily{1},raw,Dim,DD.threads.num) ;
	labBarrier
	
	%%
	[s] = initbuildRhoMean(DD);
	buildRhoMean(s,Dim);
	%%
	calcOW(d,Dim,raw,s);
	labBarrier
	%%
	metaD=d{1};
	metaD.dim=Dim;
	save metaD.mat metaD
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcOW(d,Dim,raw,MeanStuff)
	todo=d.daily.timesteps;
	T=disp_progress('init','building okubo weiss netcdfs')  ;
	for zz = 0:Dim.ws(1)-1
		T=disp_progress('show',T,Dim.ws(1),Dim.ws(1))  ;
		strt = [zz 0 0] ;
		len = [1 Dim.ws(2:3)] ;
		Tt=disp_progress('init','looping over time')  ;
		for tt = todo;
			Tt=disp_progress('show',Tt,numel(d.daily.timesteps),10)  ;
			rhoHighPass=nc_varget(d.daily.Fout{tt+1},'density',strt,len) - nc_varget(MeanStuff.Fout,'RhoMean',strt,len) ;
			%% rho gradient
			[gr.drdx,gr.drdy] = getDrhodx(rhoHighPass,raw.dx,raw.dy);
			%% velocities
			vels = getVels(raw.corio.GOverF,gr,raw.depth(zz+1));
			%% uvgrads
			uvg = UVgrads(vels,raw.dx,raw.dy);
			%% deformation
			deform = getDefo(uvg);
			%% okubo weiss
			OW = permute(okuweiss(deform),[3,1,2]);
			%             nc_varput(d.daily.OWFout{tt+1},'OkuboWeiss',OW,strt,len);
			NCcritWrite(d.daily.OWFout{tt+1},'OkuboWeiss',OW,strt,len)
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ow = okuweiss(d)
	ow = (-d.vorticity.*2+d.divergence.*2+d.stretch.*2+d.shear.*2)/2;
	ow(abs(ow)>1) = nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defo = getDefo(uvg)
	defo.vorticity = uvg.dVdx - uvg.dUdy;
	defo.divergence = uvg.dUdx + uvg.dVdy;
	defo.stretch = uvg.dUdx - uvg.dVdy;
	defo.shear = uvg.dVdx + uvg.dUdy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uvg = UVgrads(vels,dx,dy)
	%% calc U gradients
	dUdy = diff(vels.U,1,1);
	dUdx = diff(vels.U,1,2);
	dVdy = diff(vels.V,1,1);
	dVdx = diff(vels.V,1,2);
	uvg.dUdy = dUdy( [1:end, end], :)./ dy;
	uvg.dUdx = dUdx( :,[1:end, end] )./ dx;
	uvg.dVdy = dVdy( [1:end, end], :)./ dy;
	uvg.dVdx = dVdx( :,[1:end, end] )./ dx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vels = getVels(GOverF,gr,depth)
	rhoRef = 1000;
	gzOverRhoF = GOverF * depth / rhoRef;
	vels.U = -gr.drdy .* gzOverRhoF;
	vels.V = gr.drdx .* gzOverRhoF;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [drdx,drdy] = getDrhodx(rho,dx,dy)
	%% calc density gradients
	drdx = diff(rho,1,2);
	drdy = diff(rho,1,1);
	drdx = drdx(:,[1:end, end]) ./ dx;
	drdy = drdy([1:end, end],:) ./ dy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = initbuildRhoMean(DD)
	s.lims = DD.TSow.lims.inZ-1;
	s.Zsteps = s.lims(labindex,1):s.lims(labindex,2);
	s.files=DD.path.TSow.rho;
	s.Fout=[DD.path.TSow.dailyRhoName 'mean.nc'];
	initNcFile(s.Fout,'RhoMean',DD.TSow.window.size);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = initbuildRho(Dim,DD)
	dispM('building Rho',1)
	s.id = labindex;
	s.lims = DD.TSow.lims.inTime-1;
	s.timesteps = s.lims(s.id,1):s.lims(s.id,2);
	s.keys = DD.TS.keys;
	s.Ysteps=round(linspace(0,Dim.ws(2)-1,100));
	T=disp_progress('init','init Rho')  ;
	for cc =1:numel(s.Ysteps)-1
		T=disp_progress('show',T,numel(s.Ysteps)-1,100)  ;
		s.Ychunks{cc}=s.Ysteps(cc):s.Ysteps(cc+1)-1;
	end
	s.Fin = DD.path.TSow.files;
	s.dirOut=DD.path.full3d.name;
	s.Fout=DD.path.TSow.rho;
	s.OWFout=DD.path.TSow.OW;
	s.geoOut=DD.path.TSow.geo;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  buildRhoMean(s,Dim)
	oneDit = @(md) reshape(md(:),[],1);
	spmd
		rhoMean=codistributed(nan(Dim.ws(1)*Dim.ws(2)*Dim.ws(3),1));
		T=disp_progress('init','building density mean')  ;
		for ff = 1:numel(s.files)
			T=disp_progress('show',T,numel(s.files),10)  ;
			rhoMean=nansum([rhoMean,codistributed(oneDit(nc_varget(s.files{ff},'density'))),2]);
		end
		labBarrier
	end
	nc_varput(s.Fout,'rhoMean',reshape(gather(rhoMean),[Dim.ws(1),Dim.ws(2),Dim.ws(3)]),[0 0 0],[Dim.ws(1),Dim.ws(2),Dim.ws(3)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buildRho(s,raw,Dim,threads)
	oneDit = @(md) reshape(md,[],1);
	locCo=@(x) getLocalPart(codistributed(oneDit(x)));
	rdepth=raw.depth;
	rlat=raw.lat;
	ws=Dim.ws;
	spmd(threads)
		
		depth = locCo(repmat(double(rdepth),ws(2)*ws(3),1));
		lat = locCo(rlat);
		T=disp_progress('init','building density netcdfs')  ;
	end
	for tt = s.timesteps
		spmd(threads)
			T=disp_progress('show',T,numel(s.timesteps),numel(s.timesteps))  ;
			[temp,salt]=TSget(s.Fin(tt+1),s.keys,locCo);
			RHO=makeRho(salt,temp,sw_pres(depth,lat));
		end
		nc_varput(s.Fout{tt+1},'density',RHO{1},[0 0 0], [Dim.ws(1),Dim.ws(2),Dim.ws(3)]);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=makeRho(salt,temp,pres)
	rho = sw_dens(salt,temp,pres);
	rho(abs(rho>1e10) | rho==0)=nan;
	R=gop(@vertcat, rho,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,S]=TSget(FileIn,keys,locCo)
	T= locCo(nc_varget(FileIn.temp,keys.temp));
	S= locCo(nc_varget(FileIn.salt,keys.salt)*1000);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(fname,toAdd,WinSize)
	nc_create_empty(fname,'clobber');
	nc_adddim(fname,'k_index',WinSize(1));
	nc_adddim(fname,'i_index',WinSize(3));
	nc_adddim(fname,'j_index',WinSize(2));
	%% rho
	varstruct.Name = toAdd;
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'k_index','j_index','i_index' };
	nc_addvar(fname,varstruct)
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dim,raw] = setup(DD)
	ws = DD.TSow.window.size;
	FileIn = DD.path.TSow.files(1);
	keys = DD.TS.keys;
	raw.depth = nc_varget(FileIn.salt,keys.depth);
	raw.lat = nc_varget(FileIn.temp, keys.lat);
	raw.lon = nc_varget(FileIn.temp, keys.lon);
	[raw.dy,raw.dx] = getdydx( raw.lat, raw.lon);
	raw.corio = coriolisStuff(raw.lat);
	Dim.ws=ws;
	
	%% geo
	nc_varput(DD.path.TSow.geo,'depth',raw.depth);
	nc_varput(DD.path.TSow.geo,'lat',raw.lat);
	nc_varput(DD.path.TSow.geo,'lon',raw.lon);
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

