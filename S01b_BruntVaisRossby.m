%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% needs one 3D salt and temperature file each.
% integrates over depth to calculate
% -Brunt Väisälä frequency
% -Rossby Radius
% -Rossby wave first baroclinic phase speed
% NOTE: DOES NOT YET WORK FOR WINDOW SPANNING ZONAL BORDER OF DATA (140611,'yymmdd')
function S01b_BruntVaisRossby
	%% set up
	[DD,lims]=set_up;
	%% spmd
	main(DD,lims)
	%% make netcdf
	DD.path.Rossby=WriteNCfile(DD,lims);    
	%% update DD
	save_info(DD);
end
function main(DD,lims)
	if DD.debugmode
		spmd_body(DD,lims);
	else
		spmd(DD.threads.num)
			spmd_body(DD,lims);
		end
	end
end
function [DD,lims]=set_up
	%% init
	DD=initialise;
	%% set number of chunks to split large data (see input_vars.m)
	splits=DD.RossbyStuff.splits;
	%% threads
	DD.threads.num=init_threads(DD.threads.num);
	%% find temp and salt files
	[DD.path.TempSalt.salt,DD.path.TempSalt.temp]=tempsalt(DD);	
    [DD.TS.window,grids]=GetWindow(DD.path.TempSalt.salt,DD,DD.map.TS.pattern);
    [CUT]=ZonalProblem(grids,DD.TS.window);
    %% set dimension for splitting (files dont fit in memory)
	X=DD.TS.window.size.X;
	%% map chunks
	lims.data=thread_distro(splits,X) + DD.TS.window.limits.west-1;
	%% distro chunks to threads
	lims.loop=thread_distro(DD.threads.num,splits);
end
function spmd_body(DD,lims)
	id=labindex;
	%% loop over chunks
	for chnk=lims.loop(id,1):lims.loop(id,2)
		Calculations(DD,lims.data,chnk);
	end
end
function Calculations(DD,lims,chnk)
	cc=[sprintf(['%0',num2str(length(num2str(size(lims,1)))),'i'],chnk),'/',num2str(size(lims,1))];
	disp('initialising..')
	CK=initCK(DD,lims,chnk);
	%% calculate Brunt-Väisälä f and potential vorticity
	[CK.BRVA]=calcBrvaPvort(CK,cc);
	%% integrate first baroclinic rossby radius
	[CK.rossby.Ro1]=calcRossbyRadius(CK,cc);
	%% rossby wave phase speed
	[CK.rossby.c1]=calcC_one(CK,cc);
	%% save
	disp('saving..')
	saveChunk(CK,DD,chnk);
end
function nc_file_name = WriteNCfile(DD,lims)
	nc_file_name=initNC(DD);
	splits=DD.RossbyStuff.splits;
	T=disp_progress('init','creating netcdf');
	for chnk=1:splits
		T=disp_progress('disp',T,splits,100);
		%% put chunks back 2g4
		catChunks2NetCDF(DD,lims,chnk,nc_file_name)
	end
end
function	nc_file_name = initNC(DD)
	nc_file_name=[DD.path.Rossby.name, 'BVRf_all.nc'];
	overwriteornot(nc_file_name)
end
function overwriteornot(nc_file_name)
	%% test if file exists already
	try
		nc_create_empty(nc_file_name,'noclobber')
	catch me
		disp(me.message)
		reply = input('Do you want to overwrite? Y/N [Y]: ', 's');
		if isempty(reply)
			reply = 'Y';
		end
		if strcmp(reply,'Y')
			nc_create_empty(nc_file_name,'clobber')
		else
			error('exiting')
		end
	end
end
function catChunks2NetCDF(DD,lims,chnk,nc_file_name)
	%% init
	CK=loadChunk(DD,chnk);
	[Z,Y,X]=size(CK.BRVA);
	strt=lims.data(chnk,1)-DD.map.window.limits.west;
	dim.start2d  = [  0 strt];
	dim.len2d = [Y X];
	dim.start1d   = [0];
	dim.len1d = [Z];
	%% add dimensions to netcdf
	if chnk==1
		nc_adddim(nc_file_name,'depth_diff',Z);
		nc_adddim(nc_file_name,'i_index',DD.map.window.size.X	);
		nc_adddim(nc_file_name,'j_index',DD.map.window.size.Y);
	end
	%% Ro1
	varstruct.Name = 'RossbyRadius';
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'j_index','i_index' };
	if chnk==1,nc_addvar(nc_file_name,varstruct); end
	nc_varput(nc_file_name,'RossbyRadius',CK.rossby.Ro1,dim.start2d, dim.len2d);
	%% c1
	varstruct.Name = 'RossbyPhaseSpeed';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'j_index','i_index' };
	if chnk==1,nc_addvar(nc_file_name,varstruct); end
	nc_varput(nc_file_name,'RossbyPhaseSpeed',CK.rossby.c1,dim.start2d, dim.len2d);
end
function saveChunk(CK,DD,chnk) %#ok<*INUSL>
	file_out=[DD.path.TempSalt.name,'BVRf_',sprintf('%03d',chnk),'.mat'];
	save(file_out,'-struct','CK');
end
function CK=loadChunk(DD,chnk)
	file_in=[DD.path.TempSalt.name,'BVRf_',sprintf('%03d',chnk),'.mat'];
	CK=load(file_in);
end
function R=	calcRossbyRadius(CK,cc)
	disp(['integrating Rossby Radius for chunk ',cc])
	[~,YY,XX]=size(CK.BRVA);
	M.depthdiff=repmat(diff(CK.DEPTH),[1 YY XX]);
	R=abs(double((squeeze(nansum(M.depthdiff.*CK.BRVA,1))./CK.rossby.f)/pi));
end
function [c1]=calcC_one(CK,cc)
	%    c=-beta/(k^2+(1/L_r)^2) approx -beta*L^2
	disp(['applying long rossby wave disp rel for c1 for chunk ',cc])
	c1=-CK.rossby.beta.*CK.rossby.Ro1.^2;
end
function [BRVA]=calcBrvaPvort(CK,cc)
	[ZZ,YY,XX]=size(CK.TEMP);
	disp(['calculating brunt väisälä, chunk ',cc]);
	%% get full matrices for all variables
	M.depth=double(repmat(CK.DEPTH,[1,YY*XX]));
	M.lat=double(repmat(permute(CK.LAT(:),[2 1]), [ZZ,1]));
	M.pressure=double(reshape(sw_pres(M.depth(:),M.lat(:)),[ZZ,YY*XX]));
	M.salt=double(reshape(CK.SALT,[ZZ,YY*XX]));
	M.temp=double(reshape(CK.TEMP,[ZZ,YY*XX]));
	%% get brunt väisälä frequency and pot vort
	[brva,~,~]=sw_bfrq(M.salt,M.temp,M.pressure,M.lat);
	brva(brva<0)=nan;
	BRVA=sqrt(reshape(brva,[ZZ-1,YY,XX]));
end
function [CK,DD]=initCK(DD,lims,chnk)
	CK.chunk=chnk;
	dim=ncArrayDims(DD,lims,chnk);
	disp('getting temperature..')
	CK.TEMP=ChunkTemp(DD,dim);
	disp('getting salt..')
	CK.SALT=ChunkSalt(DD,dim);
	disp('getting depth..')
	CK.DEPTH=ChunkDepth(DD);
	disp('getting geo info..')
	[CK.LAT,CK.LON]=ChunkLatLon(DD,dim);
	[CK.rossby]=ChunkRossby(CK);
end
function [rossby]=ChunkRossby(CK)
	day_sid=23.9344696*60*60;
	om=2*pi/(day_sid); % frequency earth
	rossby.f=2*om*sind(CK.LAT);
	rossby.beta=2*om/earthRadius*cosd(CK.LAT);
end

function [lat,lon]=ChunkLatLon(DD,dim)
	lat=nc_varget(DD.path.TempSalt.temp,DD.map.TS.pattern.lat,dim.start1d, dim.len1d);
	lon=nc_varget(DD.path.TempSalt.temp,DD.map.TS.pattern.lon,dim.start1d, dim.len1d);
end
function depth=ChunkDepth(DD)
	depth=nc_varget(DD.path.TempSalt.salt,'depth_t');
end
function salt=ChunkSalt(DD,dim)
% 	dispNcInfo(DD.path.TempSalt.salt)
	salt=squeeze(nc_varget(DD.path.TempSalt.salt,'SALT',dim.start2d,dim.len2d));
	salt(salt==0)=nan;
	salt=salt*1000; % to salinity unit. TODO: from input vars
end
function temp=ChunkTemp(DD,dim)
% 	dispNcInfo(DD.path.TempSalt.temp)
	temp=squeeze(nc_varget(DD.path.TempSalt.temp,'TEMP',dim.start2d,dim.len2d));
	temp(temp==0)=nan;
end
function dispNcInfo(ncIn)
	%% works for the POP data...
	try
		info=nc_info(ncIn);
		disp(info.Dataset(end-1).Name);
		disp(info.Dataset(end-1).Dimension);
	catch disperr
		disp(disperr)
	end
end
function dim=ncArrayDims(DD,lims,chnk)
	j_indx_start = DD.map.window.limits.south-1;
	j_len = DD.map.window.size.Y;
	dim.start2d = [0 0 j_indx_start lims(chnk,1)-1];
	dim.len2d = 	[inf inf j_len diff(lims(chnk,:))+1];
	dim.start1d = [j_indx_start lims(chnk,1)-1];
	dim.len1d =	[j_len diff(lims(chnk,:))+1];
end
function [fileS,fileT]=tempsalt(DD)
	%% find the temp and salt files
	for kk=1:numel(DD.path.TempSalt.files)
		if ~isempty(strfind(upper(DD.path.TempSalt.files(kk).name),'SALT'))
			fileS=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name];
		end
		if ~isempty(strfind(upper(DD.path.TempSalt.files(kk).name),'TEMP'))
			fileT=[DD.path.TempSalt.name DD.path.TempSalt.files(kk).name];
		end
	end
end

