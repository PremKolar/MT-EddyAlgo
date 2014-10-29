function makePressure
	rootDir='/scratch/uni/ifmto/u300065/PUBLIC/STrhoP9495/';
	outDir = [rootDir 'pressure/'];
	[lalo,rhoFiles,sshFiles]=getData(rootDir);
	%%
	addpath(genpath('../'));
	mkdirp(outDir);
	init_threads(12);
	%%
	[Z,Y,X,dDEP,~,G]= inits(lalo);
	parfor ii=1:numel(rhoFiles)
		opDay(sshFiles,rhoFiles,ii,G,dDEP,Y,X,Z,numel(rhoFiles),outDir)
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lalo,rhoFiles,sshFiles]=getData(rD)
	rhoFiles=dir2([rD 'rho/rho_*.nc']);
	sshFiles=dir2([rD 'ssh/SSH_*.nc']);
	lld=[rD 'LatLonDepth.nc'];
	lalo.la=nc_varget(lld,'lat');
	lalo.lo=nc_varget(lld,'lon');
	lalo.de=nc_varget(lld,'depth');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Z,Y,X,dDEP,LAT,G]= inits(lalo)
	[Y,X]=  size(lalo.la);
	[Z]= numel(lalo.de) ;
	disp('dep')
	DEP=[zeros(1,Y,X); repmat(lalo.de,[1,Y,X])];
	ga=@(M) M(:);
	disp('lat')
	LAT=(permute(repmat(lalo.la,[1,1,Z]),[3,1,2]));
	disp('g')
	G=(reshape(sw_g(LAT(:),ga(DEP(2:end,:,:))),[Z,Y,X]));
	dDEP=diff(DEP,1,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opDay(sshFiles,rhoFiles,ii,G,dDEP,Y,X,Z,II,outDir)
	fprintf('%d%%\n',round(ii/II)*100);
	Frho = rhoFiles(ii).fullname;
	Srho = sshFiles(ii).fullname;
	dens  = (nc_varget(Frho,'density'));
	P_ageos     = cumsum(G.*dens.*dDEP,1);
	ssh  = repmat(permute(nc_varget(Srho,'SSH')/100,[3,1,2]),[Z,1,1]); %cm2m
	rhoZero=1000;
	P_zero         = (G .* ssh) * rhoZero;
	P = P_zero + P_ageos;
	% 	if exist(Pname,'file'), return;end
	Pname = [outDir 'pressure_' Frho(5:end)];
	ncOp(Pname,P,X,Y,Z,'pressure');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ncOp(Pname,P,X,Y,Z,fieldname)
	nc_create_empty(Pname,'noclobber');
	nc_adddim(Pname,'i_index',X);
	nc_adddim(Pname,'j_index',Y);
	nc_adddim(Pname,'k_index',Z);
	%%
	varstruct.Name = fieldname;
	varstruct.Nctype = 'single';
	varstruct.Dimension = {'k_index','j_index','i_index' };
	nc_addvar(Pname,varstruct);
	%%
	nc_varput(Pname,fieldname,single(P));
	%%
	varstruct.Name = 'depth';
	varstruct.Nctype = 'double';
	varstruct.Dimension = {'k_index'};
	nc_addvar(Pname,varstruct);
	%%
	nc_varput(Pname,'depth',dep);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S=dir2(d)
	S=dir(d);
	[base,~,~]=fileparts(S);
	for ii=1:numel(S)
		S(ii).fullname=[base S(ii).name];
	end
end
