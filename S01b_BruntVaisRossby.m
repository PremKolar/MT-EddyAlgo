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
	if DD.TS.overwrite
		main(DD,lims)
	end
	%% make netcdf
	WriteNCfile(DD,lims);
	%% update DD
	save_info(DD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,lims)
	if DD.debugmode
		spmd_body(DD,lims);
	else
		spmd(DD.threads.num)
			spmd_body(DD,lims);
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DD,lims]=set_up
	%% init
	DD=initialise;
	%% check if exists already
	[DD.path.Rossby , DD.TS.overwrite] = initNC(DD);
	%% set number of chunks to split large data (see input_vars.m)
	splits=DD.RossbyStuff.splits;
	%% threads
	DD.threads.num=init_threads(DD.threads.num);
	%% find temp and salt files
	[DD.path.TempSalt.salt,DD.path.TempSalt.temp]=tempsalt(DD);
	file=DD.path.TempSalt.salt;
	[DD.TS.window,~]=GetWindow(file,DD.map.in,DD.TS.keys);
	%% set dimension for splitting (files dont fit in memory)
	X=DD.TS.window.size.X;
	%% map chunks
	temp=thread_distro(splits,X) + DD.TS.window.limits.west-1;
	%% in case window crosses zonal bndry
	temp(temp>DD.TS.window.fullsize(2)) = temp(temp>DD.TS.window.fullsize(2)) - DD.TS.window.fullsize(2);
	%% in case one chunk crosses zonal bndry
	td=temp(:,2)-temp(:,1) < 0; % find chunk
	temp(td,1)=0; % let it start at 0
	temp(find(td)-1,2)=DD.TS.window.fullsize(2)-1; % let the one before finish at end(x-index)
	lims.data=temp;
	%% distro chunks to threads
	lims.loop=thread_distro(DD.threads.num,splits);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD,lims)
	id=labindex;
	%% loop over chunks
	for chnk=lims.loop(id,1):lims.loop(id,2)
		file_out=[DD.path.Rossby.name,'BVRf_',sprintf('%03d',chnk),'.mat'];
		if exist(file_out,'file'), disp([num2str(chnk) ' exists']);continue;end
		Calculations(DD,lims.data,chnk);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(DD)
	CK=loadChunk(DD,1);
	[dim.Z,dim.Y,dim.X]=size(CK.BRVA);
	nc_adddim(DD.path.Rossby,'depth_diff',dim.Z);
	nc_adddim(DD.path.Rossby,'i_index',DD.TS.window.size.X);
	nc_adddim(DD.path.Rossby,'j_index',DD.TS.window.size.Y);
	%% Ro1
	varstruct.Name = 'RossbyRadius';
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'j_index','i_index' };
	nc_addvar(DD.path.Rossby,varstruct)
	%% c1
	varstruct.Name = 'RossbyPhaseSpeed';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'j_index','i_index' };
	nc_addvar(DD.path.Rossby,varstruct)
	%% lat
	varstruct.Name = 'lat';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'j_index','i_index' };
	nc_addvar(DD.path.Rossby,varstruct)
	%% lon
	varstruct.Name = 'lon';
	varstruct.Nctype = 'double';
	varstruct.Dimension ={'j_index','i_index' };
	nc_addvar(DD.path.Rossby,varstruct)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CK]=initNcChunk(DD,chnk)
	CK=loadChunk(DD,chnk);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx,out]=reallocCase(DD)
	
	out=load([DD.path.cuts.name DD.path.cuts.files(1).name],'grids');
	out.lat=out.grids.lat;
	out.lon=out.grids.lon;
	out.inc.x=max(mean(diff(out.lon,1,2),1));
	out.inc.y=max(mean(diff(out.lat,1,1),2));
	[out.dim.y,out.dim.x]=size(out.lat);
	out.grids=[];
	in.lat=nc_varget(DD.path.Rossby,'lat');
	in.lon=nc_varget(DD.path.Rossby,'lon');
	
	%% get codisp'ed indeces
	lims=thread_distro(DD.threads.num,numel(in.lon));
	idx=zeros(1,numel(in.lon));
	spmd(DD.threads.num)
		idx=realloc(in,out,lims,idx);
		%% merge composite
		idxx=gop(@vertcat,idx,1);
	end
	idx=sum(idxx{1});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nc2mat(DD)
	fns=fieldnames(nc_getall(DD.path.Rossby));
	in.currVar=nc_varget(DD.path.Rossby,fns{1});
	if numel(in.currVar)~=prod(struct2array(DD.map.window.size))
		reallocIdx=true;
		[idx,~]=reallocCase(DD,DD.path.Rossby);
	end
	for fn=fns';fn=fn{1};
		%% get pop data
		in.currVar=nc_varget(nc_file_name,fn);
		out=in.currVar;
		save([DD.path.Rossby.name, fn,'.mat'],'out');
		if reallocIdx
			downscalePop(in,fn,nan(size(out.lat)),DD,idx)
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteNCfile(DD,lims)
	%% init
	splits=DD.RossbyStuff.splits;
	initNcFile(DD)
	T=disp_progress('init','creating netcdf');
	for chnk=1:splits
		T=disp_progress('disp',T,splits,100);
		%% put chunks back 2g4
		[CK,dim]=initNcChunk(DD,chnk,lims);
		catChunks2NetCDF(DD.path.Rossby,dim,CK);
	end
	%% make also mat files
	nc2mat(DD,DD.path.Rossby)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx=realloc(in,out,lims,idx)
	%% get geo info for output
	JJ=lims(labindex,1):lims(labindex,2);
	idx=getIndicesForOutMaps(in,out,JJ,idx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  downscalePop(in,fn,out,DD,idx)
	%% get geo info for output
	uni=unique(idx(idx~=0 & ~isnan(idx)));
	T=disp_progress('init','downscaling');
	system(['mv ' [DD.path.Rossby.name, fn,'.mat'] ' ' [DD.path.Rossby.name, fn,'PopSize.mat']]);
	for li=uni
		T=disp_progress('show',T,numel(uni),10);
		out(li)=nanmean(  in.currVar(idx==li));
	end
	save([DD.path.Rossby.name, fn,'.mat'],'out');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function	[outfilename,overwrite] = initNC(DD)
	outfilename=[DD.path.Rossby.name, 'BVRf_all.nc'];
	overwrite=NCoverwriteornot(outfilename);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function catChunks2NetCDF(DD,dim,CK)
	nc_varput(DD.path.Rossby,'RossbyRadius',CK.rossby.Ro1,dim.start2d, dim.len2d);
	nc_varput(DD.path.Rossby,'RossbyPhaseSpeed',CK.rossby.c1,dim.start2d, dim.len2d);
	nc_varput(DD.path.Rossby,'lat',CK.lat	,dim.start2d, dim.len2d);
	nc_varput(DD.path.Rossby,'lon',CK.lon	,dim.start2d, dim.len2d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveChunk(CK,DD,chnk) %#ok<*INUSL>
	file_out=[DD.path.Rossby.name,'BVRf_',sprintf('%03d',chnk),'.mat'];
	save(file_out,'-struct','CK');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CK=loadChunk(DD,chnk)
	file_in=[DD.path.Rossby.name,'BVRf_',sprintf('%03d',chnk),'.mat'];
	CK=load(file_in);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=	calcRossbyRadius(CK,cc)
	disp(['integrating Rossby Radius for chunk ',cc])
	[~,YY,XX]=size(CK.BRVA);
	M.depthdiff=repmat(diff(CK.DEPTH),[1 YY XX]);
	R=abs(double((squeeze(nansum(M.depthdiff.*CK.BRVA,1))./CK.rossby.f)/pi));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c1]=calcC_one(CK,cc)
	%    c=-beta/(k^2+(1/L_r)^2) approx -beta*L^2
	disp(['applying long rossby wave disp rel for c1 for chunk ',cc])
	c1=-CK.rossby.beta.*CK.rossby.Ro1.^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BRVA]=calcBrvaPvort(CK,cc)
	[ZZ,YY,XX]=size(CK.TEMP);
	disp(['calculating brunt väisälä, chunk ',cc]);
	%% get full matrices for all variables
	M.depth=double(repmat(CK.DEPTH,[1,YY*XX]));
	M.lat=double(repmat(permute(CK.lat(:),[2 1]), [ZZ,1]));
	M.pressure=double(reshape(sw_pres(M.depth(:),M.lat(:)),[ZZ,YY*XX]));
	M.salt=double(reshape(CK.SALT,[ZZ,YY*XX]));
	M.temp=double(reshape(CK.TEMP,[ZZ,YY*XX]));
	%% get brunt väisälä frequency and pot vort
	[brva,~,~]=sw_bfrq(M.salt,M.temp,M.pressure,M.lat);
	brva(brva<0)=nan;
	BRVA=sqrt(reshape(brva,[ZZ-1,YY,XX]));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CK,DD]=initCK(DD,lims,chunk)
	CK.chunk=chunk;
	CK.dim=ncArrayDims(DD,lims,chunk);
	disp('getting temperature..')
	CK.TEMP=ChunkTemp(DD,CK.dim);
	disp('getting salt..')
	CK.SALT=ChunkSalt(DD,CK.dim);
	disp('getting depth..')
	CK.DEPTH=ChunkDepth(DD);
	disp('getting geo info..')
	[CK.lat,CK.lon]=ChunkLatLon(DD,CK.dim);
	[CK.rossby]=ChunkRossby(CK);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rossby]=ChunkRossby(CK)
	day_sid=23.9344696*60*60;
	om=2*pi/(day_sid); % frequency earth
	rossby.f=2*om*sind(CK.lat);
	rossby.beta=2*om/earthRadius*cosd(CK.lat);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lat,lon]=ChunkLatLon(DD,dim)
	lat=nc_varget(DD.path.TempSalt.temp,DD.TS.keys.lat,dim.start1d, dim.len1d);
	lon=nc_varget(DD.path.TempSalt.temp,DD.TS.keys.lon,dim.start1d, dim.len1d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function depth=ChunkDepth(DD)
	depth=nc_varget(DD.path.TempSalt.salt,'depth_t');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function salt=ChunkSalt(DD,dim)
	salt=squeeze(nc_varget(DD.path.TempSalt.salt,'SALT',dim.start2d,dim.len2d));
	salt(salt==0)=nan;
	salt=salt*1000; % to salinity unit. TODO: from input vars
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function temp=ChunkTemp(DD,dim)
	temp=squeeze(nc_varget(DD.path.TempSalt.temp,'TEMP',dim.start2d,dim.len2d));
	temp(temp==0)=nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dim=ncArrayDims(DD,lims,chnk)
	j_indx_start = DD.TS.window.limits.south-1;
	j_len = DD.TS.window.size.Y;
	dim.start2d = [0 0 j_indx_start lims(chnk,1)-1];
	dim.len2d = 	[inf inf j_len diff(lims(chnk,:))+1];
	dim.start1d = [j_indx_start lims(chnk,1)-1];
	dim.len1d =	[j_len diff(lims(chnk,:))+1];
	%% new indeces for output nc file
	xlens=diff(squeeze(lims(:,2)));
	xlens(xlens<0)= xlens(xlens<0) + DD.TS.window.fullsize(2);
	newxstart=sum(xlens(1:chnk-1))+1 -1;
	dim.new.start =[0 newxstart];
	dim.new.len =  dim.len1d;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

