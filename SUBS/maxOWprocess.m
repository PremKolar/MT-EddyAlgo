%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 17-Jul-2014 23:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  maxOWprocess
	dbstop if error
	load DD
	load metaData
	xdghn
	NC=initNC(metaData,DD.path.Rossby.name);
	spmdBcalc(NC);
	
	%
	%
	%     %% focus on strong neg. okubo weiss
	%     clean=OWall.ow < 5*nanmean(OWall.ow(:));
	%     OWall.ow(~clean) = nan;
	%     OWall.zi(~clean) = nan;
	%     %% ow weighted mean of zi
	%     OWall.owSum      = repmat(nansum(OWall.ow,1),[NC.S.T,1,1]);
	%     OWall.ziWeighted = OWall.ow.*OWall.zi./OWall.owSum;
	%     OWall.meaned.z           = squeeze(nansum(OWall.ziWeighted, 1));
	%     OWall.meaned.ow          = squeeze(nanmean(OWall.ow, 1));
	%     flgd=~squeeze(nansum(clean,1));
	%     [y,x]=size(OWall.meaned.z);
	%     [Xq,Yq]=meshgrid(1:x,1:y);
	%     Xfl=Xq;Yfl=Yq;
	%     Xfl(flgd)=[];
	%     Yfl(flgd)=[];
	%     vq = griddata(Xfl,Yfl,OWall.meaned.z(~flgd),Xq,Yq);
	%     OWall.ziIntrl=round(smooth2a(NeighbourValue(isnan(vq),vq),10));
	%     %    pcolor(vqn);
	%     %    colorbar;
	%
	%     allOW=OWall.ow;
	%     depthOW=OWall.depth;
	%     ziOW=OWall.zi;
	%     ziIntrp=OWall.ziIntrl;
	%     ziWeighted=OWall.ziWeighted;
	%
	%
	%     save('allOW.mat','allOW','-v7.3')
	%     save('zi.mat','ziOW','-v7.3')
	%     save('ziItnrp.mat','ziIntrp','-v7.3')
	%     save('ziWeighted.mat','ziWeighted','-v7.3')
	%
	%
	%
	%
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NC=initNC(metaD,outdir)
	NC.geo=metaD.geoOut;
	NC.files(numel(metaD.OWFout)).n=struct;
	NC.outdir=outdir;
	finishedF=cell(size(metaD.OWFout));
	for cc=1:numel(metaD.OWFout)
		[dr,fi,ex]=fileparts(metaD.OWFout{cc}) ;
		finishedF{cc}=[dr '/' fi ex ];
	end
	[NC.files(:).n] = deal(finishedF{:});
	NC.timeStep = datenum('0815','mmdd');
	NC.S=cell2struct(num2cell(getfield(nc_getvarinfo(NC.files(1).n,'OkuboWeiss'),'Size')),{'Z','Y','X'},2);
	NC.S.T = numel(NC.files);
	NC.Sname=cell2struct(getfield(nc_getvarinfo(NC.files(1).n,'OkuboWeiss'),'Dimension'),{'Z','Y','X'},2);
	
	%%
	S=cell2mat(struct2cell(NC.S))';
	NC.new.dimName = {'t_index','j_index','i_index' };
	NC.new.dimNum  = S([4 2 3]);
	NC.new.minOW.varName		 =		'minOW';
	NC.new.minOWzi.varName	 =		'zOfMinOW';
	NC.new.minOW.fileName	 =		[outdir 'minOW.nc'];
	NC.new.minOWzi.fileName  =		[outdir 'zOfminOW.nc'];
	%% init
	NC.iniNewNC = @(n,f,D,Dn) initNcFile(n.(f).fileName,n.(f).varName,D,Dn);
	NC.iniNewNC(NC.new,'minOWzi',NC.new.dimNum,NC.new.dimName);
	NC.iniNewNC(NC.new,'minOW',  NC.new.dimNum,NC.new.dimName);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function daily=initDaily(NC,tt)
	daily.minOWzi.varName = NC.new.minOWzi.varName;
	daily.minOW.varName   = NC.new.minOW.varName;
	daily.minOWzi.fileName = [NC.outdir sprintf('zOfminOW_%04d.nc',tt)];
	daily.minOW.fileName   = [dirbase sprintf('minOW_%04d.nc',tt)];
	NC.iniNewNC(daily,'minOWzi',NC.new.dimNum(2:end),NC.new.dimName(2:end));
	NC.iniNewNC(daily,'minOW',  NC.new.dimNum(2:end),NC.new.dimName(2:end));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmdBcalc(NC)
	f=funcs;
	ncPut=@(n,f,data)  nc_varput(n.(f).fileName ,n.(f).varName,data);
	ncPutBig=@(n,f,data,t)  nc_varput(n.(f).fileName ,n.(f).varName,data,[t,0,0],[1 inf inf]);
	%% get bathymetry
	NCf1=NC.files(1).n;
	[~,~,~,bath]=getBathym(nc_varget(NCf1,'OkuboWeiss'));
	%%
	for tt=1:NC.S.T
		daily=initDaily(NC,tt);
		%% get min in z
		log10data = f.log10OW(nc_varget(NC.files(tt).n,'OkuboWeiss'));
		[owMin,MinZi]=spmdBlck(log10data,bath,f);
		%% write daily
		ncPut(daily,'minOWzi',MinZi);
		ncPut(daily,'minOW',owMin);
		%% put to big files too
		ncPutBig(NC.new,'minOWzi',MinZi,tt-1);
		ncPutBig(NC.new,'minOW',owMin,tt-1);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z,Y,X,bath]=getBathym(OW)
	[Z,Y,X]=size(OW);
	OW2d=reshape(OW,[Z,Y*X]);
	[~,bathUpdown]=min(isnan(flipud(OW2d)),[],1);
	bath=reshape( Z-bathUpdown + 1, [Y,X]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=funcs
	f.ncvOne = @(A,dim) getLocalPart(codistributed(A,codistributor1d(dim)));
	f.nanminFrom2toFloor = @(OW,bath) nanmin(OW(2:bath-1,:,:),[], 1);
	f.gCat = @(a,dim) gcat(squeeze(a),dim,1);
	f.log10OW = @(OW,dummy) log10(-prep4log(OW,dummy));
end
function [OW]=prep4log(OW,dummy)
	OW(isnan(OW) || isinf(OW) || OW>=0)=dummy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [owMin,MinZi]=spmdBlck(data,bath,f)
	spmd
		%% min in z
		mydata=f.ncvOne(data,3);
		mybath=f.ncvOne(bath,2);
		[owMin,MinZi]=f.nanminFrom2toFloor(mydata,mybath);
		MinZi=f.gCat(MinZi-1,2); % correct for (2: ...)
		owMin=f.gCat(owMin,2);
	end
	MinZi=MinZi{1};
	owMin=owMin{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initNcFile(fname,toAdd,WinSize,dimName)
	nc_create_empty(fname,'clobber');
	varstruct.Name = toAdd;
	varstruct.Nctype = 'double';
	for ww=1:numel(WinSize)
		nc_adddim(fname,dimName{ww},WinSize(ww));
	end
	varstruct.Dimension = dimName;
	nc_addvar(fname,varstruct)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








