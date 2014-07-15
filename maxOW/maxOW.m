%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOW
	%% set up
	[DD]=maxOWsetUp;
	%% spmd
	main(DD)
	%% make netcdf
	WriteNCfile(DD);
	%% update DD
	save_info(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD)
	if DD.debugmode
		spmd_body(DD);
	else
		spmd(DD.threads.num)
			spmd_body(DD);
			disp_progress('conclude');
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD)
	id=labindex;
	lims=DD.RossbyStuff.lims;
	%%
	[CKpre]=preInitCK(DD);
	%% loop over chunks
	for ff=0:numel(DD.path.TSow)-1
		for chnk=lims.loop(id,1):lims.loop(id,2)
			Calculations(DD,chnk,ff,CKpre);
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Calculations(DD,chnk,ff,CKpre)
	%% init
	[CKnow,cc]=init(DD,chnk,ff);
	%% merge
	CK=mergeStruct2(CKpre,CKnow);
	clear CKpre CKnow;
	%% calculate Brunt-Väisälä f and potential vorticity
	[CK.pres]=calcPres(CK,cc);
	%% OW
	[CK.OW]=calcOW(CK,cc);
	CK.pres=[];
	%% save
	disp('saving..')
	saveChunk(CK,DD.path.Rossby.name,chnk,ff);
	%----------------------------------------------------------------------
	function [CK,cc]=init(DD,chnk,ff)
		lims=DD.RossbyStuff.lims.data;
		cc=[sprintf(['%0',num2str(length(num2str(size(lims,1)))),'i'],chnk),'/',num2str(size(lims,1))];
		disp('initialising..')
		CK=initCK(DD,chnk,ff);
	end
	%----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OW =	calcOW(CK,cc)
	dispM(['getting okubo weiss',cc],1)
	%% p gradient
	[gr.dpdx,gr.dpdy]=getDpdx(CK.pres,CK.DX,CK.DY);
	%% velocities
	vels=getVels(CK.corio,gr);
	clear gr;
	%% uvgrads
	uvg=UVgrads(vels,DX,DY);
	clear vels;
	%% deformation
	deform=getDefo(uvg);
	clear uvg;
	%% okubo weiss
	OW=okuweiss(deform);
	%% ---------------------------------------------------------------------
	function OW=okuweiss(d)
		OW=(-d.vorticity.*2+d.divergence.*2+d.stretch.*2+d.shear.*2)/2;
	end
	%-----------------------------------------------------------------------
	function defo=getDefo(uvg)
		defo.vorticity = uvg.dVdx - uvg.dUdy;
		defo.divergence= uvg.dUdx + uvg.dVdy;
		defo.stretch   = uvg.dUdx - uvg.dVdy;
		defo.shear     = uvg.dVdx + uvg.dUdy;
	end
	%-----------------------------------------------------------------------
	function uvg=UVgrads(vels,DX,DY)
		%% calc U gradients
		dUdy=diff(vels.U,1,2);
		dUdx=diff(vels.U,1,3);
		dVdy=diff(vels.V,1,2);
		dVdx=diff(vels.V,1,3);
		[Z,~,~]=size(vels.U);
		uvg.dUdy= dUdy(:, [1:end, end], :)  ./ vertstack(DY,Z);
		uvg.dUdx= dUdx(:, :,[1:end, end] )  ./ vertstack(DX,Z);
		uvg.dVdy= dVdy(:, [1:end, end], :)  ./ vertstack(DY,Z);
		uvg.dVdx= dVdx(:, :,[1:end, end] )  ./ vertstack(DX,Z);
	end
	%-----------------------------------------------------------------------
	function vels=getVels(cor,gr)
		GOF=vertstack(cor.GOverF,Z);
		vels.U=-GOF.*gr.dpdy;
		vels.V= GOF.*gr.dpdx;
		clear GOF;
		vels.absUV=hypot(abs(vels.U),abs(vels.V));
	end
	%-----------------------------------------------------------------------
	function [dpdx,dpdy]=getDpdx(p,DX,DY)
		%% calc ssh gradients
		[Z,~,~]=size(p);
		dpdx=diff(p,1,3);
		dpdy=diff(p,1,2);
		dpdx=dpdx(:,:,[1:end, end]) ./ vertstack(DX,Z);
		dpdy=dpdy(:,[1:end, end],:) ./ vertstack(DY,Z);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pressure]=calcPres(CK,cc)
	[ZZ,YY,XX]=size(CK.TEMP);
	dispM(['calculating pressure, chunk ',cc]);
	%% get full matrices for all variables
	M.depth=double(repmat(CK.depth,[1,YY*XX]));
	M.lat=double(repmat(permute(CK.lat(:),[2 1]), [ZZ,1]));
	pressure=double(reshape(sw_pres(M.depth(:),M.lat(:)),[ZZ,YY,XX]));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CK]=initCK(DD,chunk,ff)
	CK.chunk=chunk;
	disp('getting temperature..')
	CK.TEMP=ChunkTemp(DD,CK.dim,ff+1);
	disp('getting salt..')
	CK.SALT=ChunkSalt(DD,CK.dim+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CK]=preInitCK(DD)
	CK.dim=ncArrayDims(DD,chunk,1);
	disp('getting depth..')
	CK.depth=ChunkDepth(DD);
	disp('getting geo info..')
	[CK.lat,CK.lon]=ChunkLatLon(DD,CK.dim);
	[CK.DY,CK.DX]=ChunkDYDX(CK.lat,CK.lon);
	disp('getting coriolis stuff..')
	[CK.rossby]=ChunkRossby(CK);
	CK.corio=coriolisStuff(CK.lat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=coriolisStuff(lat)
	%% omega
	out.Omega=angularFreqEarth;
	%% f
	out.f=2*out.Omega*sind(lat);
	%% beta
	out.beta=2*out.Omega/earthRadius*cosd(lat);
	%% gravity
	out.g=sw_g(lat,zeros(size(lat)));
	%% g/f
	out.GOverF=out.g./out.f;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rossby]=ChunkRossby(CK)
	day_sid=23.9344696*60*60;
	om=2*pi/(day_sid); % frequency earth
	rossby.f=2*om*sind(CK.lat);
	rossby.beta=2*om/earthRadius*cosd(CK.lat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DY,DX]=ChunkDYDX(lat,lon)
	%% grid increment sizes
	DY=deg2rad(abs(diff(double(lat),1,1)))*earthRadius;
	DX=deg2rad(abs(diff(double(lon),1,2)))*earthRadius.*cosd(lat(:,1:end-1));
	%% append one line/row to have identical size as other fields
	DY=DY([1:end end],:);
	DX=DX(:,[1:end end]);
	%% correct 360° crossings
	seamcrossflag=DX>100*median(DX(:));
	DX(seamcrossflag)=abs(DX(seamcrossflag) - 2*pi*earthRadius.*cosd(lat(seamcrossflag)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveChunk(CK,RossbyDir,cc,ff)  %#ok<INUSL>
	file_out=[RossbyDir,'OW_',sprintf('%03d',ff),'_',sprintf('%03d',cc),'.mat'];
	save(file_out,'-struct','CK');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CK=loadChunk(RossbyDir,chnk,ff)
	file_in=[RossbyDir,'OW_',sprintf('%03d',ff),'_',sprintf('%03d',chnk),'.mat'];
	CK=load(file_in);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dim=ncArrayDims(DD,chnk,ff)
	lims=DD.RossbyStuff.lims.data;
	j_indx_start = DD.TS.window.limits.south-1;
	j_len = DD.TS.window.size.Y;
	dim.start2d = [ff 0 j_indx_start lims(chnk,1)-1];
	dim.len2d = 	[1 inf j_len diff(lims(chnk,:))+1];
	dim.start1d = [j_indx_start lims(chnk,1)-1];
	dim.len1d =	[j_len diff(lims(chnk,:))+1];
	%% new indeces for output nc file
	xlens=diff(squeeze(lims(:,2)));
	xlens(xlens<0)= xlens(xlens<0) + DD.TS.window.fullsize(2);
	newxstart=sum(xlens(1:chnk-1))+1 -1;
	dim.new.start.fourD =[ff 0 0 newxstart];
	dim.new.len.fourD =  dim.len2d;
	dim.new.start.z =[0];
	dim.new.len.z =  [inf];
	dim.new.start.twoD =[0 newxstart];
	dim.new.len.twoD =  dim.len1d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(DD)
	CK=loadChunk(DD.path.Rossby.name,1,1);
	[dim.Z,dim.Y,dim.X]=size(CK.pres);
	dim.t=numel(DD.path.TSow);
	nc_adddim(DD.path.Rossby.NCfile,'k_index',dim.Z);
	nc_adddim(DD.path.Rossby.NCfile,'i_index',DD.TS.window.size.X);
	nc_adddim(DD.path.Rossby.NCfile,'j_index',DD.TS.window.size.Y);
	nc_adddim(DD.path.Rossby.NCfile,'t_index',dim.t);
	%% OW
	varstruct.Name = 'Okubo-Weiss';
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'t_index','k_index','j_index','i_index' };
	nc_addvar(DD.path.Rossby.NCfile,varstruct)
	%% depth
	varstruct.Name = 'depth';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'k_index'};
	nc_addvar(DD.path.Rossby.NCfile,varstruct)
	%% lat
	varstruct.Name = 'lat';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'j_index','i_index' };
	nc_addvar(DD.path.Rossby.NCfile,varstruct)
	%% lon
	varstruct.Name = 'lon';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'j_index','i_index' };
	nc_addvar(DD.path.Rossby.NCfile,varstruct)
	%% time
	varstruct.Name = 'time';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'t_index'};
	nc_addvar(DD.path.Rossby.NCfile,varstruct)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function catChunks2NetCDF(file,CK,ff)
	start= CK.dim.new.start;
	len  = CK.dim.new.len;
	nc_varput(file,'Okubo-Weiss',CK.OW, start.fourD,    len.fourD);
	nc_varput(file,'depth',CK.depth,        start.z,        len.z);
	nc_varput(file,'lat', CK.lat,			    start.twoD,  len.twoD);
	nc_varput(file,'lon', CK.lon,           start.twoD,  len.twoD);
	nc_varput(file,'time', now + ff,        [ff],             [1]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteNCfile(DD)
	%% init
	splits=DD.parameters.RossbySplits;
	initNcFile(DD);
	Tf=disp_progress('init','looping over time-steps');
	fnum=numel(DD.path.TSow);
	for ff=0:fnum-1
		Tf=disp_progress('disp',Tf,fnum,100);
		%% cat chunks
		T=disp_progress('init','creating netcdf');
		for chnk=1:splits
			T=disp_progress('disp',T,splits,6);
			%% put chunks back 2g4
			CK=loadChunk(DD.path.Rossby.name,chnk,ff);
			catChunks2NetCDF(DD.path.Rossby.NCfile,CK,ff);
		end
		%% make also mat files
		nc2mat(DD);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lat,lon]=ChunkLatLon(DD,dim)
	lat=nc_varget(DD.path.TSow(1).temp,DD.TS.keys.lat,dim.start1d, dim.len1d);
	lon=nc_varget(DD.path.TSow(1).temp,DD.TS.keys.lon,dim.start1d, dim.len1d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function depth=ChunkDepth(DD)
	depth=nc_varget(DD.path.TSow(1).salt,'depth_t');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function salt=ChunkSalt(DD,dim,ff)
	salt=squeeze(nc_varget(DD.path.TSow(ff).salt,'TEMP',dim.start2d,dim.len2d));
	salt(salt==0)=nan;
	salt=salt*1000; % to salinity unit. TODO: from input vars
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp=ChunkTemp(DD,dim,ff)
	temp=squeeze(nc_varget(DD.path.TSow(ff).temp,'TEMP',dim.start2d,dim.len2d));
	temp(temp==0)=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc2matSave(DD,fn)
	%% get pop data
	out=nc_varget(DD.path.Rossby.NCfile,fn); %#ok<NASGU>
	save([DD.path.Rossby.name, fn,'.mat'],'out');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc2mat(DD)
	%% test for remap
	fns=ncfieldnames(DD.path.Rossby.NCfile)';
	%% save 2 mats
	for fn=fns;
		nc2matSave(DD,fn{1})
	end
end
